"""Convergence analysis for grid independence study."""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
import math


@dataclass
class LevelResult:
    """Results for a single mesh level."""
    level_name: str
    num_cells: int
    metric_value: float
    metric_min: Optional[float] = None
    metric_max: Optional[float] = None
    run_time: Optional[float] = None  # seconds


@dataclass
class ConvergenceResult:
    """Result of convergence comparison between two levels."""
    coarse_level: str
    fine_level: str
    coarse_value: float
    fine_value: float
    percent_change: float
    converged: bool
    threshold: float


@dataclass
class StudyAnalysis:
    """Complete analysis of a grid study."""
    level_results: List[LevelResult]
    convergence_results: List[ConvergenceResult]
    is_converged: bool
    converged_level: Optional[str] = None
    recommended_level: Optional[str] = None
    gci_fine: Optional[float] = None  # Grid Convergence Index
    extrapolated_value: Optional[float] = None  # Richardson extrapolation
    stop_reason: Optional[str] = None  # Reason for stopping (adaptive mode)


class GridAnalyzer:
    """Analyzes grid independence study results."""

    def __init__(self, threshold: float = 1.0):
        """
        Args:
            threshold: Convergence threshold in percent
        """
        self.threshold = threshold

    def check_convergence(
        self, coarse: LevelResult, fine: LevelResult
    ) -> tuple[float, bool]:
        """
        Quick check if two consecutive levels have converged.

        Args:
            coarse: Coarser mesh level result
            fine: Finer mesh level result

        Returns:
            (percent_change, converged)
        """
        if fine.metric_value == 0:
            percent_change = float("inf") if coarse.metric_value != 0 else 0.0
        else:
            percent_change = (
                abs(coarse.metric_value - fine.metric_value)
                / abs(fine.metric_value)
                * 100
            )
        return percent_change, percent_change < self.threshold

    def analyze(self, level_results: List[LevelResult]) -> StudyAnalysis:
        """
        Perform full grid independence analysis.

        Args:
            level_results: List of results ordered from coarse to fine

        Returns:
            Complete analysis with convergence status
        """
        if len(level_results) < 2:
            raise ValueError("Need at least 2 mesh levels for convergence analysis")

        # Calculate convergence between consecutive levels
        convergence_results = []
        for i in range(len(level_results) - 1):
            coarse = level_results[i]
            fine = level_results[i + 1]

            conv = self._calculate_convergence(coarse, fine)
            convergence_results.append(conv)

        # Determine overall convergence
        is_converged = False
        converged_level = None

        for conv in convergence_results:
            if conv.converged:
                is_converged = True
                converged_level = conv.coarse_level  # First level that achieves convergence
                break

        # Calculate GCI if we have at least 3 levels
        gci_fine = None
        extrapolated_value = None

        if len(level_results) >= 3:
            try:
                gci_fine, extrapolated_value = self._calculate_gci(level_results[-3:])
            except Exception:
                pass  # GCI calculation can fail for non-monotonic convergence

        # Determine recommended level
        recommended_level = self._determine_recommended_level(
            level_results, convergence_results, is_converged
        )

        return StudyAnalysis(
            level_results=level_results,
            convergence_results=convergence_results,
            is_converged=is_converged,
            converged_level=converged_level,
            recommended_level=recommended_level,
            gci_fine=gci_fine,
            extrapolated_value=extrapolated_value,
        )

    def _calculate_convergence(
        self, coarse: LevelResult, fine: LevelResult
    ) -> ConvergenceResult:
        """Calculate convergence between two mesh levels."""
        # percent_change = |coarse - fine| / fine * 100
        # Using fine as reference (assumed closer to true value)
        if fine.metric_value == 0:
            percent_change = float("inf") if coarse.metric_value != 0 else 0
        else:
            percent_change = (
                abs(coarse.metric_value - fine.metric_value)
                / abs(fine.metric_value)
                * 100
            )

        converged = percent_change < self.threshold

        return ConvergenceResult(
            coarse_level=coarse.level_name,
            fine_level=fine.level_name,
            coarse_value=coarse.metric_value,
            fine_value=fine.metric_value,
            percent_change=percent_change,
            converged=converged,
            threshold=self.threshold,
        )

    def _calculate_gci(
        self, levels: List[LevelResult]
    ) -> tuple[float, float]:
        """
        Calculate Grid Convergence Index using Richardson extrapolation.

        Based on Roache's method (1994).

        Args:
            levels: Three consecutive levels [coarse, medium, fine]

        Returns:
            (GCI_fine, extrapolated_value)
        """
        if len(levels) != 3:
            raise ValueError("GCI calculation requires exactly 3 levels")

        f1 = levels[0].metric_value  # coarsest
        f2 = levels[1].metric_value  # medium
        f3 = levels[2].metric_value  # finest

        n1 = levels[0].num_cells
        n2 = levels[1].num_cells
        n3 = levels[2].num_cells

        # Representative grid spacing (assuming 3D, proportional to N^(-1/3))
        h1 = n1 ** (-1 / 3)
        h2 = n2 ** (-1 / 3)
        h3 = n3 ** (-1 / 3)

        # Refinement ratios
        r21 = h1 / h2
        r32 = h2 / h3

        # Apparent order of convergence
        epsilon32 = f3 - f2
        epsilon21 = f2 - f1

        if epsilon32 == 0 or epsilon21 == 0:
            raise ValueError("Cannot calculate order - zero difference between levels")

        # Iterative solution for p (apparent order)
        # Using fixed-point iteration
        s = math.copysign(1, epsilon32 / epsilon21)

        if s * epsilon32 / epsilon21 <= 0:
            raise ValueError("Non-monotonic convergence - cannot calculate GCI")

        p = abs(math.log(abs(epsilon32 / epsilon21)) / math.log(r21))

        # Richardson extrapolated value
        f_ext = f3 + (f3 - f2) / (r32**p - 1)

        # GCI for fine grid (Fs = 1.25 as recommended safety factor)
        Fs = 1.25
        e_a = abs((f3 - f2) / f3)  # approximate relative error
        gci_fine = Fs * e_a / (r32**p - 1) * 100  # in percent

        return gci_fine, f_ext

    def _determine_recommended_level(
        self,
        level_results: List[LevelResult],
        convergence_results: List[ConvergenceResult],
        is_converged: bool,
    ) -> str:
        """Determine recommended mesh level based on analysis."""
        if not is_converged:
            # Not converged - recommend finest available
            return level_results[-1].level_name

        # Find first convergent level (balance of accuracy and cost)
        for i, conv in enumerate(convergence_results):
            if conv.converged:
                # Recommend the coarser level of the convergent pair
                return conv.coarse_level

        return level_results[-1].level_name


def format_summary_table(analysis: StudyAnalysis, metric_name: str = "T_avg") -> str:
    """
    Format a compact summary table showing grid resolution, metric, and change %.

    Example output:
    ┌─────────────────┬────────────┬────────────┬──────────┬────────┐
    │ Level           │ Cells      │ T_avg [K]  │ Δ [%]    │ Status │
    ├─────────────────┼────────────┼────────────┼──────────┼────────┤
    │ L1_coarse       │     50,000 │   330.1234 │      -   │   -    │
    │ L2_medium       │    150,000 │   332.5678 │   0.74   │  PASS  │
    │ L3_fine         │    400,000 │   332.8901 │   0.10   │  PASS  │
    └─────────────────┴────────────┴────────────┴──────────┴────────┘
    """
    # Build change % lookup
    change_map = {}
    for cr in analysis.convergence_results:
        change_map[cr.fine_level] = (cr.percent_change, cr.converged)

    # Column widths
    w_level = 17
    w_cells = 12
    w_metric = 12
    w_change = 10
    w_status = 8

    lines = []

    # Top border
    lines.append(f"┌{'─'*w_level}┬{'─'*w_cells}┬{'─'*w_metric}┬{'─'*w_change}┬{'─'*w_status}┐")

    # Header
    lines.append(
        f"│ {'Level':<{w_level-2}} │ {'Cells':>{w_cells-2}} │ {metric_name:>{w_metric-2}} │ {'Δ [%]':>{w_change-2}} │ {'Status':^{w_status-2}} │"
    )

    # Header separator
    lines.append(f"├{'─'*w_level}┼{'─'*w_cells}┼{'─'*w_metric}┼{'─'*w_change}┼{'─'*w_status}┤")

    # Data rows
    for lr in analysis.level_results:
        cells_str = f"{lr.num_cells:,}"
        metric_str = f"{lr.metric_value:.4f}"

        if lr.level_name in change_map:
            pct, converged = change_map[lr.level_name]
            change_str = f"{pct:.2f}"
            status_str = "PASS" if converged else "FAIL"
        else:
            change_str = "-"
            status_str = "-"

        lines.append(
            f"│ {lr.level_name:<{w_level-2}} │ {cells_str:>{w_cells-2}} │ {metric_str:>{w_metric-2}} │ {change_str:>{w_change-2}} │ {status_str:^{w_status-2}} │"
        )

    # Bottom border
    lines.append(f"└{'─'*w_level}┴{'─'*w_cells}┴{'─'*w_metric}┴{'─'*w_change}┴{'─'*w_status}┘")

    # Summary line
    status = "✓ CONVERGED" if analysis.is_converged else "✗ NOT CONVERGED"
    lines.append(f"\nResult: {status} (threshold: {analysis.convergence_results[0].threshold}%)")

    if analysis.recommended_level:
        lines.append(f"Recommended: {analysis.recommended_level}")

    return "\n".join(lines)


def format_analysis_table(analysis: StudyAnalysis) -> str:
    """Format analysis results as a detailed text table."""
    lines = []

    # Header
    lines.append("=" * 70)
    lines.append("GRID INDEPENDENCE STUDY RESULTS")
    lines.append("=" * 70)
    lines.append("")

    # Level results table
    lines.append("Mesh Level Results:")
    lines.append("-" * 70)
    lines.append(f"{'Level':<15} {'Cells':>12} {'Value':>15} {'Min':>12} {'Max':>12}")
    lines.append("-" * 70)

    for lr in analysis.level_results:
        min_str = f"{lr.metric_min:.4f}" if lr.metric_min else "N/A"
        max_str = f"{lr.metric_max:.4f}" if lr.metric_max else "N/A"
        lines.append(
            f"{lr.level_name:<15} {lr.num_cells:>12,} {lr.metric_value:>15.4f} "
            f"{min_str:>12} {max_str:>12}"
        )

    lines.append("")

    # Convergence results table
    lines.append("Convergence Analysis:")
    lines.append("-" * 70)
    lines.append(
        f"{'Comparison':<25} {'Coarse':>12} {'Fine':>12} {'Change %':>10} {'Status':>10}"
    )
    lines.append("-" * 70)

    for cr in analysis.convergence_results:
        comparison = f"{cr.coarse_level} → {cr.fine_level}"
        status = "PASS" if cr.converged else "FAIL"
        lines.append(
            f"{comparison:<25} {cr.coarse_value:>12.4f} {cr.fine_value:>12.4f} "
            f"{cr.percent_change:>10.2f} {status:>10}"
        )

    lines.append("")

    # Summary
    lines.append("Summary:")
    lines.append("-" * 70)
    lines.append(f"Convergence threshold: {analysis.convergence_results[0].threshold}%")
    lines.append(f"Grid independent: {'YES' if analysis.is_converged else 'NO'}")

    if analysis.converged_level:
        lines.append(f"First converged level: {analysis.converged_level}")

    if analysis.recommended_level:
        lines.append(f"Recommended level: {analysis.recommended_level}")

    if analysis.gci_fine is not None:
        lines.append(f"GCI (fine): {analysis.gci_fine:.2f}%")

    if analysis.extrapolated_value is not None:
        lines.append(f"Richardson extrapolated value: {analysis.extrapolated_value:.4f}")

    lines.append("=" * 70)

    return "\n".join(lines)
