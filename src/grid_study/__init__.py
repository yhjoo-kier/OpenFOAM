# Grid Independence Study Framework for OpenFOAM CHT Simulations
"""
Framework for automated grid independence studies.

Usage (Python API):
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

Usage (CLI):
    python -m grid_study run cases/heatsink_water_cht_steady -o output_dir
    python -m grid_study init -o config.json
    python -m grid_study run-config config.json

Convenience function:
    from grid_study import run_grid_study

    analysis = run_grid_study(
        base_case="cases/heatsink_water_cht_steady",
        output_dir="grid_study_output",
        metric_patch="heat_source",
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
__version__ = "0.1.0"
