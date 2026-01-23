#!/usr/bin/env python3
"""Command-line interface for grid study framework."""

import argparse
import sys
from pathlib import Path


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Grid Independence Study Framework for OpenFOAM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings (4 mesh levels)
  python -m grid_study run cases/heatsink_water_cht_steady

  # Run with custom output directory and threshold
  python -m grid_study run cases/my_case -o my_study -t 0.5

  # Generate config file for customization
  python -m grid_study init -o my_study_config.json

  # Run from config file
  python -m grid_study run-config my_study_config.json
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # Run command
    run_parser = subparsers.add_parser("run", help="Run grid study")
    run_parser.add_argument("base_case", type=str, help="Path to base OpenFOAM case")
    run_parser.add_argument(
        "-o", "--output", type=str, default="grid_study_output", help="Output directory"
    )
    run_parser.add_argument(
        "-n", "--name", type=str, default="grid_study", help="Study name"
    )
    run_parser.add_argument(
        "-p", "--patch", type=str, default="heat_source", help="Patch for metric extraction"
    )
    run_parser.add_argument(
        "-f", "--field", type=str, default="T", help="Field to extract"
    )
    run_parser.add_argument(
        "-r", "--region", type=str, default="solid", help="Region (for multi-region)"
    )
    run_parser.add_argument(
        "-t", "--threshold", type=float, default=1.0, help="Convergence threshold (%%)"
    )
    run_parser.add_argument(
        "--levels",
        type=int,
        default=4,
        choices=[2, 3, 4, 5],
        help="Number of mesh levels",
    )

    # Init command - generate config file
    init_parser = subparsers.add_parser("init", help="Generate config file template")
    init_parser.add_argument(
        "-o", "--output", type=str, default="grid_study_config.json", help="Output file"
    )
    init_parser.add_argument(
        "-c", "--case", type=str, help="Base case path to include in config"
    )

    # Run from config command
    config_parser = subparsers.add_parser("run-config", help="Run from config file")
    config_parser.add_argument("config_file", type=str, help="Path to config JSON file")

    # Analyze command - analyze existing results
    analyze_parser = subparsers.add_parser(
        "analyze", help="Analyze existing study results"
    )
    analyze_parser.add_argument("study_dir", type=str, help="Path to study output directory")

    args = parser.parse_args()

    if args.command == "run":
        cmd_run(args)
    elif args.command == "init":
        cmd_init(args)
    elif args.command == "run-config":
        cmd_run_config(args)
    elif args.command == "analyze":
        cmd_analyze(args)
    else:
        parser.print_help()


def cmd_run(args):
    """Run grid study command."""
    from .config import GridStudyConfig, MeshLevel
    from .study import GridStudy

    base_case = Path(args.base_case)
    if not base_case.exists():
        print(f"Error: Base case not found: {base_case}")
        sys.exit(1)

    # Get subset of default mesh levels based on --levels argument
    all_levels = GridStudyConfig.get_default_mesh_levels()
    if args.levels < len(all_levels):
        # Select levels evenly distributed from coarse to fine
        indices = [
            int(i * (len(all_levels) - 1) / (args.levels - 1))
            for i in range(args.levels)
        ]
        mesh_levels = [all_levels[i] for i in indices]
    else:
        mesh_levels = all_levels

    config = GridStudyConfig(
        base_case_path=base_case,
        study_name=args.name,
        output_dir=Path(args.output),
        mesh_levels=mesh_levels,
        metric_patch=args.patch,
        metric_field=args.field,
        metric_region=args.region,
        convergence_threshold=args.threshold,
    )

    study = GridStudy(config)
    analysis = study.run()
    reports = study.generate_reports(analysis)

    print("\n" + "=" * 60)
    print("STUDY COMPLETE")
    print("=" * 60)
    print(f"Reports saved to: {config.output_dir / 'results'}")


def cmd_init(args):
    """Generate config file template."""
    from .config import GridStudyConfig

    config = GridStudyConfig(
        base_case_path=Path(args.case) if args.case else Path("cases/your_case"),
        study_name="grid_study",
        output_dir=Path("grid_study_output"),
    )

    output_path = Path(args.output)
    config.save(output_path)
    print(f"Config template saved to: {output_path}")
    print("Edit this file to customize mesh levels and settings.")


def cmd_run_config(args):
    """Run grid study from config file."""
    from .config import GridStudyConfig
    from .study import GridStudy

    config_path = Path(args.config_file)
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    config = GridStudyConfig.load(config_path)

    if not config.base_case_path.exists():
        print(f"Error: Base case not found: {config.base_case_path}")
        sys.exit(1)

    study = GridStudy(config)
    analysis = study.run()
    reports = study.generate_reports(analysis)

    print("\n" + "=" * 60)
    print("STUDY COMPLETE")
    print("=" * 60)
    print(f"Reports saved to: {config.output_dir / 'results'}")


def cmd_analyze(args):
    """Analyze existing study results."""
    from .config import GridStudyConfig
    from .analyzer import GridAnalyzer, LevelResult
    from .reporter import StudyReporter
    import json

    study_dir = Path(args.study_dir)
    config_path = study_dir / "config.json"
    results_json = study_dir / "results" / "*_report.json"

    # Find the JSON report
    json_reports = list((study_dir / "results").glob("*_report.json"))
    if not json_reports:
        print(f"Error: No JSON report found in {study_dir / 'results'}")
        sys.exit(1)

    # Load the report
    with open(json_reports[0]) as f:
        data = json.load(f)

    # Reconstruct level results
    level_results = [
        LevelResult(
            level_name=lr["level_name"],
            num_cells=lr["num_cells"],
            metric_value=lr["metric_value"],
            metric_min=lr.get("metric_min"),
            metric_max=lr.get("metric_max"),
            run_time=lr.get("run_time"),
        )
        for lr in data["level_results"]
    ]

    # Load config
    config = GridStudyConfig.load(config_path) if config_path.exists() else None
    threshold = data["config"].get("threshold", 1.0)

    # Re-analyze
    analyzer = GridAnalyzer(threshold)
    analysis = analyzer.analyze(level_results)

    # Print results
    from .analyzer import format_analysis_table

    print(format_analysis_table(analysis))


if __name__ == "__main__":
    main()
