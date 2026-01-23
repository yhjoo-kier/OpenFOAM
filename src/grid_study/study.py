"""Main grid study orchestration."""

import time
from pathlib import Path
from typing import List, Optional, Callable
from dataclasses import dataclass

from .config import GridStudyConfig, MeshLevel
from .case_manager import CaseManager
from .mesh_generator import generate_mesh
from .runner import SimulationRunner, run_mesh_conversion, RunResult
from .extractor import ResultExtractor, MetricResult
from .analyzer import GridAnalyzer, LevelResult, StudyAnalysis, format_summary_table
from .reporter import StudyReporter


@dataclass
class StudyProgress:
    """Progress information for callbacks."""
    current_level: int
    total_levels: int
    level_name: str
    stage: str  # "mesh", "convert", "run", "extract"
    message: str


class GridStudy:
    """
    Main orchestrator for grid independence studies.

    Example usage:
        config = GridStudyConfig(
            base_case_path=Path("cases/heatsink_water_cht_steady"),
            study_name="heatsink_grid_study",
            output_dir=Path("grid_study_output"),
        )

        study = GridStudy(config)
        analysis = study.run()
        study.generate_reports(analysis)
    """

    def __init__(
        self,
        config: GridStudyConfig,
        progress_callback: Optional[Callable[[StudyProgress], None]] = None,
    ):
        """
        Initialize grid study.

        Args:
            config: Study configuration
            progress_callback: Optional callback for progress updates
        """
        self.config = config
        self.progress_callback = progress_callback or self._default_progress
        self.case_manager = CaseManager(config)
        self.runner = SimulationRunner(config.solver_app)
        self.analyzer = GridAnalyzer(config.convergence_threshold)
        self.reporter = StudyReporter(config.output_dir)

        self.level_results: List[LevelResult] = []
        self.run_results: List[RunResult] = []

    def _default_progress(self, progress: StudyProgress) -> None:
        """Default progress callback - prints to console."""
        print(
            f"[{progress.current_level}/{progress.total_levels}] "
            f"{progress.level_name}: {progress.stage} - {progress.message}"
        )

    def _report_progress(
        self,
        current: int,
        total: int,
        level_name: str,
        stage: str,
        message: str,
    ) -> None:
        """Report progress through callback."""
        self.progress_callback(
            StudyProgress(
                current_level=current,
                total_levels=total,
                level_name=level_name,
                stage=stage,
                message=message,
            )
        )

    def run(self) -> StudyAnalysis:
        """
        Run the complete grid independence study.

        Returns:
            StudyAnalysis with all results
        """
        print("=" * 60)
        print("GRID INDEPENDENCE STUDY")
        print(f"Study: {self.config.study_name}")
        print(f"Base case: {self.config.base_case_path}")
        print(f"Levels: {len(self.config.mesh_levels)}")
        print("=" * 60)
        print()

        # Setup directories
        self.case_manager.setup_study_directory()

        # Save configuration
        self.config.save(self.config.output_dir / "config.json")

        total_levels = len(self.config.mesh_levels)

        for i, mesh_level in enumerate(self.config.mesh_levels, 1):
            print(f"\n{'='*60}")
            print(f"LEVEL {i}/{total_levels}: {mesh_level.name}")
            print(f"{'='*60}")

            try:
                level_result = self._run_level(i, total_levels, mesh_level)
                self.level_results.append(level_result)

                print(f"  Result: {level_result.metric_value:.4f}")
                print(f"  Cells: {level_result.num_cells:,}")

            except Exception as e:
                print(f"  ERROR: {e}")
                raise

        # Analyze results
        print(f"\n{'='*60}")
        print("ANALYZING RESULTS")
        print(f"{'='*60}")

        analysis = self.analyzer.analyze(self.level_results)

        # Print summary
        self._print_summary(analysis)

        return analysis

    def _run_level(
        self,
        current: int,
        total: int,
        mesh_level: MeshLevel,
    ) -> LevelResult:
        """Run a single mesh level."""
        level_start_time = time.time()

        # 1. Create case directory
        self._report_progress(current, total, mesh_level.name, "setup", "Creating case")
        case_dir = self.case_manager.create_case(mesh_level)

        # Update controlDict with configured iterations
        self.case_manager.update_control_dict(case_dir, self.config.num_iterations)

        # 2. Generate mesh
        self._report_progress(current, total, mesh_level.name, "mesh", "Generating mesh")
        msh_file = self.config.output_dir / "meshes" / f"{mesh_level.name}.msh"
        msh_file.parent.mkdir(parents=True, exist_ok=True)

        mesh_info = generate_mesh(msh_file, mesh_level)
        print(f"  Mesh generated: {mesh_info['num_cells']:,} cells")

        # 3. Convert mesh to OpenFOAM format
        self._report_progress(
            current, total, mesh_level.name, "convert", "Converting to OpenFOAM format"
        )
        log_dir = self.config.output_dir / "logs"

        success = run_mesh_conversion(case_dir, msh_file, log_dir)
        if not success:
            raise RuntimeError(f"Mesh conversion failed for {mesh_level.name}")

        # 4. Run simulation
        self._report_progress(
            current, total, mesh_level.name, "run", "Running simulation"
        )
        run_result = self.runner.run(case_dir, log_dir, mesh_level.name)
        self.run_results.append(run_result)

        if not run_result.success:
            raise RuntimeError(
                f"Simulation failed for {mesh_level.name}: {run_result.error_message}"
            )

        print(f"  Simulation completed in {run_result.elapsed_time:.1f}s")

        # 5. Extract metrics
        self._report_progress(
            current, total, mesh_level.name, "extract", "Extracting results"
        )
        extractor = ResultExtractor(case_dir)

        metric_result = extractor.get_patch_average(
            field=self.config.metric_field,
            patch=self.config.metric_patch,
            region=self.config.metric_region,
        )

        mesh_info_full = extractor.get_mesh_info()
        total_cells = mesh_info_full.get("total_cells", mesh_info["num_cells"])

        level_elapsed = time.time() - level_start_time

        return LevelResult(
            level_name=mesh_level.name,
            num_cells=total_cells,
            metric_value=metric_result.value,
            metric_min=metric_result.min_value,
            metric_max=metric_result.max_value,
            run_time=level_elapsed,
        )

    def _print_summary(self, analysis: StudyAnalysis) -> None:
        """Print analysis summary to console using compact table format."""
        print()
        # Use metric name with unit for table header
        metric_label = f"{self.config.metric_name} [K]"
        print(format_summary_table(analysis, metric_label))

        if analysis.extrapolated_value:
            print(f"Richardson extrapolated: {analysis.extrapolated_value:.4f} K")

    def generate_reports(self, analysis: StudyAnalysis) -> dict:
        """
        Generate all report formats.

        Args:
            analysis: Completed study analysis

        Returns:
            Dictionary with paths to generated reports
        """
        config_summary = {
            "base_case": str(self.config.base_case_path),
            "metric_name": self.config.metric_name,
            "metric_patch": self.config.metric_patch,
            "metric_field": self.config.metric_field,
            "metric_region": self.config.metric_region,
            "threshold": self.config.convergence_threshold,
        }

        reports = {}

        # Text report
        reports["text"] = self.reporter.generate_text_report(
            analysis, self.config.study_name, config_summary
        )
        print(f"Text report: {reports['text']}")

        # JSON report
        reports["json"] = self.reporter.generate_json_report(
            analysis, self.config.study_name, config_summary
        )
        print(f"JSON report: {reports['json']}")

        # CSV report
        reports["csv"] = self.reporter.generate_csv_report(
            analysis, self.config.study_name
        )
        print(f"CSV report: {reports['csv']}")

        # Plot
        plot_path = self.reporter.generate_plot(
            analysis, self.config.study_name, self.config.metric_name
        )
        if plot_path:
            reports["plot"] = plot_path
            print(f"Plot: {reports['plot']}")

        return reports


def run_grid_study(
    base_case: str,
    output_dir: str = "grid_study_output",
    study_name: str = "grid_study",
    metric_patch: str = "heat_source",
    metric_field: str = "T",
    metric_region: str = "solid",
    threshold: float = 1.0,
) -> StudyAnalysis:
    """
    Convenience function to run a grid study.

    Args:
        base_case: Path to base OpenFOAM case
        output_dir: Output directory for results
        study_name: Name for the study
        metric_patch: Patch to extract metric from
        metric_field: Field to extract (e.g., "T", "p")
        metric_region: Region for multi-region cases
        threshold: Convergence threshold in percent

    Returns:
        StudyAnalysis with results
    """
    config = GridStudyConfig(
        base_case_path=Path(base_case),
        study_name=study_name,
        output_dir=Path(output_dir),
        metric_patch=metric_patch,
        metric_field=metric_field,
        metric_region=metric_region,
        convergence_threshold=threshold,
    )

    study = GridStudy(config)
    analysis = study.run()
    study.generate_reports(analysis)

    return analysis
