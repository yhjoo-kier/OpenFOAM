"""Report generation for grid study results."""

import json
from pathlib import Path
from datetime import datetime
from typing import Optional

from .analyzer import StudyAnalysis, LevelResult


class StudyReporter:
    """Generates reports for grid independence studies."""

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.results_dir = output_dir / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)

    def generate_text_report(
        self,
        analysis: StudyAnalysis,
        study_name: str,
        config_summary: dict,
    ) -> Path:
        """Generate a text report."""
        report_path = self.results_dir / f"{study_name}_report.txt"

        lines = []
        lines.append("=" * 80)
        lines.append(f"GRID INDEPENDENCE STUDY REPORT")
        lines.append(f"Study: {study_name}")
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("=" * 80)
        lines.append("")

        # Configuration summary
        lines.append("CONFIGURATION")
        lines.append("-" * 80)
        lines.append(f"Base case: {config_summary.get('base_case', 'N/A')}")
        lines.append(f"Metric: {config_summary.get('metric_name', 'N/A')}")
        lines.append(f"Patch: {config_summary.get('metric_patch', 'N/A')}")
        lines.append(f"Region: {config_summary.get('metric_region', 'N/A')}")
        lines.append(f"Convergence threshold: {config_summary.get('threshold', 1.0)}%")
        lines.append("")

        # Mesh levels
        lines.append("MESH LEVELS")
        lines.append("-" * 80)
        lines.append(f"{'Level':<20} {'Cells':>15} {'Run Time':>15}")
        lines.append("-" * 80)

        for lr in analysis.level_results:
            time_str = f"{lr.run_time:.1f}s" if lr.run_time else "N/A"
            lines.append(f"{lr.level_name:<20} {lr.num_cells:>15,} {time_str:>15}")

        lines.append("")

        # Results
        lines.append("RESULTS")
        lines.append("-" * 80)
        lines.append(f"{'Level':<20} {'Value':>15} {'Min':>15} {'Max':>15}")
        lines.append("-" * 80)

        for lr in analysis.level_results:
            min_str = f"{lr.metric_min:.4f}" if lr.metric_min else "N/A"
            max_str = f"{lr.metric_max:.4f}" if lr.metric_max else "N/A"
            lines.append(
                f"{lr.level_name:<20} {lr.metric_value:>15.4f} {min_str:>15} {max_str:>15}"
            )

        lines.append("")

        # Convergence analysis
        lines.append("CONVERGENCE ANALYSIS")
        lines.append("-" * 80)

        for cr in analysis.convergence_results:
            status = "✓ CONVERGED" if cr.converged else "✗ NOT CONVERGED"
            lines.append(f"{cr.coarse_level} → {cr.fine_level}:")
            lines.append(f"  Change: {cr.percent_change:.4f}% (threshold: {cr.threshold}%)")
            lines.append(f"  Status: {status}")
            lines.append("")

        # Summary
        lines.append("SUMMARY")
        lines.append("-" * 80)
        lines.append(f"Grid independent: {'YES' if analysis.is_converged else 'NO'}")

        if analysis.converged_level:
            lines.append(f"First converged level: {analysis.converged_level}")

        if analysis.recommended_level:
            lines.append(f"Recommended mesh level: {analysis.recommended_level}")

        if analysis.gci_fine is not None:
            lines.append(f"Grid Convergence Index (fine): {analysis.gci_fine:.4f}%")

        if analysis.extrapolated_value is not None:
            lines.append(
                f"Richardson extrapolated value: {analysis.extrapolated_value:.4f}"
            )

        lines.append("")
        lines.append("=" * 80)

        report_path.write_text("\n".join(lines))
        return report_path

    def generate_json_report(
        self,
        analysis: StudyAnalysis,
        study_name: str,
        config_summary: dict,
    ) -> Path:
        """Generate a JSON report for programmatic access."""
        report_path = self.results_dir / f"{study_name}_report.json"

        data = {
            "study_name": study_name,
            "timestamp": datetime.now().isoformat(),
            "config": config_summary,
            "level_results": [
                {
                    "level_name": lr.level_name,
                    "num_cells": lr.num_cells,
                    "metric_value": lr.metric_value,
                    "metric_min": lr.metric_min,
                    "metric_max": lr.metric_max,
                    "run_time": lr.run_time,
                }
                for lr in analysis.level_results
            ],
            "convergence_results": [
                {
                    "coarse_level": cr.coarse_level,
                    "fine_level": cr.fine_level,
                    "coarse_value": cr.coarse_value,
                    "fine_value": cr.fine_value,
                    "percent_change": cr.percent_change,
                    "converged": cr.converged,
                    "threshold": cr.threshold,
                }
                for cr in analysis.convergence_results
            ],
            "summary": {
                "is_converged": analysis.is_converged,
                "converged_level": analysis.converged_level,
                "recommended_level": analysis.recommended_level,
                "gci_fine": analysis.gci_fine,
                "extrapolated_value": analysis.extrapolated_value,
            },
        }

        with open(report_path, "w") as f:
            json.dump(data, f, indent=2)

        return report_path

    def generate_csv_report(
        self,
        analysis: StudyAnalysis,
        study_name: str,
    ) -> Path:
        """Generate a CSV report for spreadsheet analysis."""
        report_path = self.results_dir / f"{study_name}_results.csv"

        lines = [
            "level_name,num_cells,metric_value,metric_min,metric_max,run_time_s,percent_change,converged"
        ]

        for i, lr in enumerate(analysis.level_results):
            # Get convergence info if available (compares to next finer level)
            pct_change = ""
            converged = ""

            # Find convergence result where this level is the coarse one
            for cr in analysis.convergence_results:
                if cr.coarse_level == lr.level_name:
                    pct_change = f"{cr.percent_change:.4f}"
                    converged = "yes" if cr.converged else "no"
                    break

            min_val = lr.metric_min if lr.metric_min else ""
            max_val = lr.metric_max if lr.metric_max else ""
            run_time = lr.run_time if lr.run_time else ""

            lines.append(
                f"{lr.level_name},{lr.num_cells},{lr.metric_value:.6f},"
                f"{min_val},{max_val},{run_time},{pct_change},{converged}"
            )

        report_path.write_text("\n".join(lines))
        return report_path

    def generate_plot(
        self,
        analysis: StudyAnalysis,
        study_name: str,
        metric_name: str = "Metric Value",
    ) -> Optional[Path]:
        """Generate convergence plot."""
        try:
            import matplotlib.pyplot as plt
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
        except ImportError:
            print("matplotlib not available, skipping plot generation")
            return None

        plot_path = self.results_dir / f"{study_name}_convergence.png"

        # Prepare data
        levels = [lr.level_name for lr in analysis.level_results]
        cells = [lr.num_cells for lr in analysis.level_results]
        values = [lr.metric_value for lr in analysis.level_results]

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Plot 1: Metric vs number of cells
        ax1.plot(cells, values, "bo-", linewidth=2, markersize=8)
        ax1.set_xlabel("Number of Cells")
        ax1.set_ylabel(metric_name)
        ax1.set_title("Metric vs Mesh Size")
        ax1.set_xscale("log")
        ax1.grid(True, alpha=0.3)

        # Add level labels
        for i, (c, v, l) in enumerate(zip(cells, values, levels)):
            ax1.annotate(
                l.replace("_", "\n"),
                (c, v),
                textcoords="offset points",
                xytext=(0, 10),
                ha="center",
                fontsize=8,
            )

        # Add extrapolated value if available
        if analysis.extrapolated_value:
            ax1.axhline(
                y=analysis.extrapolated_value,
                color="r",
                linestyle="--",
                label=f"Extrapolated: {analysis.extrapolated_value:.4f}",
            )
            ax1.legend()

        # Plot 2: Percent change vs level transition
        if analysis.convergence_results:
            transitions = [
                f"{cr.coarse_level}\n→\n{cr.fine_level}"
                for cr in analysis.convergence_results
            ]
            changes = [cr.percent_change for cr in analysis.convergence_results]
            colors = [
                "green" if cr.converged else "red"
                for cr in analysis.convergence_results
            ]

            bars = ax2.bar(transitions, changes, color=colors, alpha=0.7)
            ax2.axhline(
                y=analysis.convergence_results[0].threshold,
                color="orange",
                linestyle="--",
                label=f"Threshold ({analysis.convergence_results[0].threshold}%)",
            )
            ax2.set_xlabel("Level Transition")
            ax2.set_ylabel("Percent Change (%)")
            ax2.set_title("Convergence Between Levels")
            ax2.legend()
            ax2.grid(True, alpha=0.3, axis="y")

            # Add value labels on bars
            for bar, val in zip(bars, changes):
                height = bar.get_height()
                ax2.annotate(
                    f"{val:.2f}%",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        plt.close()

        return plot_path
