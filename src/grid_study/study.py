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

        In standard mode, runs all predefined mesh levels.
        In adaptive mode, continues refining until convergence or limits reached.

        Returns:
            StudyAnalysis with all results
        """
        mode_str = "ADAPTIVE" if self.config.adaptive_mode else "STANDARD"
        print("=" * 60)
        print(f"GRID INDEPENDENCE STUDY ({mode_str} MODE)")
        print(f"Study: {self.config.study_name}")
        print(f"Base case: {self.config.base_case_path}")
        if self.config.adaptive_mode:
            print(f"Convergence threshold: {self.config.convergence_threshold}%")
            print(f"Max cells: {self.config.max_cells:,}")
            print(f"Max levels: {self.config.max_levels}")
        else:
            print(f"Levels: {len(self.config.mesh_levels)}")
        print("=" * 60)
        print()

        # Setup directories
        self.case_manager.setup_study_directory()

        # Save configuration
        self.config.save(self.config.output_dir / "config.json")

        if self.config.adaptive_mode:
            analysis = self._run_adaptive()
        else:
            analysis = self._run_standard()

        # Print summary
        self._print_summary(analysis)

        return analysis

    def _run_standard(self) -> StudyAnalysis:
        """Run standard mode with predefined mesh levels."""
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

        return self.analyzer.analyze(self.level_results)

    def _run_adaptive(self) -> StudyAnalysis:
        """Run adaptive mode, refining until convergence or limits reached."""
        # Start with predefined levels if available, otherwise use defaults
        mesh_levels = list(self.config.mesh_levels) if self.config.mesh_levels else []

        if not mesh_levels:
            # Create initial coarse level
            mesh_levels = [
                MeshLevel(
                    name="L1_initial",
                    mesh_factor=2.0,
                    bl_first_height=0.0005,
                    bl_growth_ratio=1.2,
                    bl_num_layers=0,
                )
            ]

        stop_reason = None
        level_idx = 0

        while True:
            level_idx += 1

            # Check max levels limit
            if level_idx > self.config.max_levels:
                stop_reason = f"max_levels_reached ({self.config.max_levels})"
                print(f"\n⚠ Stopping: Maximum levels reached ({self.config.max_levels})")
                break

            # Get or generate current mesh level
            if level_idx <= len(mesh_levels):
                mesh_level = mesh_levels[level_idx - 1]
            else:
                # Generate next level based on previous
                prev_level = mesh_levels[-1]
                mesh_level = self.config.generate_next_level(prev_level, level_idx)
                mesh_levels.append(mesh_level)

            print(f"\n{'='*60}")
            print(f"LEVEL {level_idx}: {mesh_level.name} (mesh_factor={mesh_level.mesh_factor:.3f})")
            print(f"{'='*60}")

            try:
                level_result = self._run_level(level_idx, "?", mesh_level)
                self.level_results.append(level_result)

                print(f"  Result: {level_result.metric_value:.4f}")
                print(f"  Cells: {level_result.num_cells:,}")

                # Check max cells limit for next level
                if level_result.num_cells > self.config.max_cells:
                    stop_reason = f"max_cells_exceeded ({level_result.num_cells:,} > {self.config.max_cells:,})"
                    print(f"\n⚠ Stopping: Maximum cells exceeded ({level_result.num_cells:,})")
                    break

                # Check convergence (need at least 2 levels)
                if len(self.level_results) >= 2:
                    prev_result = self.level_results[-2]
                    curr_result = self.level_results[-1]
                    pct_change, converged = self.analyzer.check_convergence(prev_result, curr_result)

                    status = "✓ CONVERGED" if converged else f"Δ={pct_change:.2f}%"
                    print(f"  Convergence: {status}")

                    if converged:
                        stop_reason = f"converged (Δ={pct_change:.2f}% < {self.config.convergence_threshold}%)"
                        print(f"\n✓ Convergence achieved at {mesh_level.name}!")
                        break

                # Check max runtime per level
                if self.config.max_runtime_per_level:
                    if level_result.run_time and level_result.run_time > self.config.max_runtime_per_level:
                        stop_reason = f"max_runtime_exceeded ({level_result.run_time:.0f}s > {self.config.max_runtime_per_level}s)"
                        print(f"\n⚠ Stopping: Level runtime exceeded limit ({level_result.run_time:.0f}s)")
                        break

            except Exception as e:
                print(f"  ERROR: {e}")
                stop_reason = f"error: {str(e)}"
                break

        # Analyze results
        print(f"\n{'='*60}")
        print("ANALYZING RESULTS")
        print(f"{'='*60}")

        if len(self.level_results) < 2:
            raise RuntimeError("Adaptive mode requires at least 2 levels to complete")

        analysis = self.analyzer.analyze(self.level_results)
        analysis.stop_reason = stop_reason

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

        if analysis.stop_reason:
            print(f"Stop reason: {analysis.stop_reason}")

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
    adaptive: bool = False,
    max_cells: int = 2_000_000,
    max_levels: int = 10,
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
        adaptive: Enable adaptive refinement mode
        max_cells: Maximum cells limit (adaptive mode)
        max_levels: Maximum levels limit (adaptive mode)

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
        adaptive_mode=adaptive,
        max_cells=max_cells,
        max_levels=max_levels,
    )

    study = GridStudy(config)
    analysis = study.run()
    study.generate_reports(analysis)

    return analysis
