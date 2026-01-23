# Grid Independence Study Framework for OpenFOAM CHT Simulations
"""
Framework for automated grid independence studies.

Usage (Python API - Standard Mode):
    from grid_study import GridStudy, GridStudyConfig
    from pathlib import Path

    config = GridStudyConfig(
        base_case_path=Path("cases/heatsink_water_cht_steady"),
        study_name="heatsink_grid_study",
        output_dir=Path("grid_study_output"),
    )

    study = GridStudy(config)
    analysis = study.run()
    reports = study.generate_reports(analysis)

Usage (Python API - Adaptive Mode):
    config = GridStudyConfig(
        base_case_path=Path("cases/heatsink_water_cht_steady"),
        study_name="adaptive_study",
        adaptive_mode=True,           # Enable adaptive refinement
        convergence_threshold=1.0,    # Stop when Δ < 1%
        max_cells=2_000_000,          # Stop if cells exceed this
        max_levels=10,                # Stop after this many levels
        refinement_ratio=0.7,         # Each level: mesh_factor *= 0.7
    )

    study = GridStudy(config)
    analysis = study.run()  # Continues until converged or limits reached

Convenience function:
    from grid_study import run_grid_study

    # Standard mode
    analysis = run_grid_study(
        base_case="cases/heatsink_water_cht_steady",
        metric_patch="heat_source",
    )

    # Adaptive mode
    analysis = run_grid_study(
        base_case="cases/heatsink_water_cht_steady",
        adaptive=True,
        max_cells=1_000_000,
    )
"""

from .config import GridStudyConfig, MeshLevel
from .study import GridStudy, run_grid_study
from .analyzer import GridAnalyzer, StudyAnalysis, LevelResult, format_summary_table

__all__ = [
    "GridStudyConfig",
    "MeshLevel",
    "GridStudy",
    "run_grid_study",
    "GridAnalyzer",
    "StudyAnalysis",
    "LevelResult",
    "format_summary_table",
]
__version__ = "0.2.0"
