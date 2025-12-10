"""
Staggered circular finned-tube bank performance model
======================================================

This module implements correlation-based performance calculations for a
staggered circular tube – circular fin heat exchanger under a fixed
gas-side pumping power (or pressure-drop) constraint.

The implementation follows:
  - **Heat transfer**: Briggs & Young (1963) correlation
  - **Pressure drop**: ESDU 86022 (1986) correlation

참조 문서:
  - 고정 펌핑 동력 제약 조건 하에서의 엇갈림 배열 원형 핀-관 열교환기 최적화에 관한 심층 연구 보고서
  - 스태거드 원형 관-원형 핀 열교환기의 최적 배열 연구
  - docs/ESDU_86022_상관식_정리.md

The user-facing variable names are mapped to the symbols in the original
correlations as follows (all lengths in [m]):

Geometry of the tube bank external volume
-----------------------------------------
  - W_bank : Tube bank width (frontal width normal to flow)
  - H_bank : Tube bank height (tube axis direction, fin stack direction)
  - L_bank : Tube bank length in flow direction (depth, number of tube rows)

Tube & fin geometry
-------------------
  - d_o    : Tube outer diameter          (Briggs/Robinson: d_o)
  - S_T    : Transverse tube pitch        (Robinson: S_t)
             Pitch between tube center-lines in the spanwise direction
             (perpendicular to flow, within the frontal plane).
  - S_L    : Longitudinal tube pitch      (Robinson: S_l)
             Pitch between tube center-lines in the flow direction.
  - S_fin  : Clear fin spacing s          (Briggs: s)
             여기서는 **두 휜 사이의 비어 있는 간격(clear spacing)** 으로 정의한다.
             Fin pitch (center-to-center) = t_fin + S_fin.
  - t_fin  : Fin thickness                (Briggs: t)
  - h_fin  : Radial fin height            (Briggs: h_f)
             h_fin = (D_fin - d_o) / 2, where D_fin is the fin outer diameter.

Fluid properties (gas side, default: air at 0°C)
------------------------------------------------
  - rho_gas : Density          [kg/m^3]
  - cp_gas  : Specific heat    [J/kg-K]
  - mu_gas  : Dynamic viscosity [Pa·s]
  - k_gas   : Thermal conductivity [W/m-K]
  - Pr_gas  : Prandtl number   (= mu_gas * cp_gas / k_gas)

Dimensionless groups
--------------------
  - Re_d  : Reynolds number based on tube outer diameter d_o
            Re_d = rho_gas * V_max * d_o / mu_gas
            where V_max is the maximum velocity in the minimum flow area.

  - Nu    : Nusselt number based on d_o
            Nu = h * d_o / k_gas

  - j     : Colburn j-factor
            j = Nu / (Re_d * Pr_gas**(1/3))

  - f     : Friction factor K_f (ESDU 86022 definition, per tube row)

Briggs & Young heat-transfer correlation (high-fin circular tubes)
------------------------------------------------------------------
  Nu = 0.134 * Re_d**0.681 * Pr_gas**(1/3) * (s / h_f)**0.2 * (s / t)**0.1134

with
  s   = S_fin (clear fin spacing)  [m]
  h_f = h_fin                      [m]
  t   = t_fin                      [m]

ESDU 86022 friction correlation (high-fin staggered tube bank)
---------------------------------------------------------------
ESDU 86022 (1986)는 high-fin staggered tube bank에 대한 산업 표준 상관식으로,
Robinson & Briggs (1966)보다 더 넓은 핀 간격 범위를 다룬다.

압력강하:
  ΔP = (K_acc + N_rows * K_f) * (1/2) * ρ * V_max²

마찰계수:
  K_f = 4.567 * Re_d^(-0.242) * A_increase^0.504
        * (S_T/d_o)^(-0.376) * (S_L/d_o)^(-0.546)

가속손실:
  K_acc = 1 + σ²

with
  A_increase = A_total / A_bare_tube  (핀 포함 표면적 / 베어튜브 표면적)
  σ = A_min / A_fr  (flow area contraction ratio)

적용 범위 (ESDU 86022):
  참조: docs/ESDU 86022 압력강하 상관식 적용 범위 - Gemini.md

  - Re_d (V_max 기준): 10² ~ 10⁵ (층류 말단 ~ 완전 난류)
  - S_T/d_o: 1.1 ~ 4.0 (가장 높은 신뢰도 범위)
  - S_L/d_o: 1.1 ~ 3.0 (가장 높은 신뢰도 범위)
  - Fin density: 4 ~ 11 fpi (s ≈ 2.0 ~ 6.1 mm for t_fin = 0.25 mm)
  - d_o: 9.5 ~ 50.8 mm
  - h_fin: 8.5 ~ 15.9 mm (본 코드에서는 7.0 mm까지 완화 적용)
  - d_f/d_o: 1.2 ~ 2.4

참고: docs/ESDU 86022 압력강하 상관식 적용 범위 - Gemini.md

Fixed pumping power constraint
------------------------------
For a given pumping power target P_pump_target [W], we assume

  P_pump = (m_dot / rho_gas) * ΔP

where ΔP is the total pressure drop across the tube bank on the gas side.

We solve for the frontal velocity V_fr (or equivalently m_dot) such that

  P_pump(V_fr) = P_pump_target

using a 1D root-finding (bisection) method.

Area definitions (for reporting)
--------------------------------
We report three surface areas:

  - A_fin_total  : Total fin external area (both sides + perimeter) [m^2]
  - A_tube_total : Total external tube area, approximated as
                   π * d_o * H_bank * N_tubes (ignoring shading by fins).
                   이 값은 핀에 의해 가려지는 부분을 포함하므로,
                   \"실제 공기와 접촉하는 면적\"보다 클 수 있다.
  - A_total      : A_fin_total + A_tube_total

The local convective coefficient from the Briggs & Young correlation is

  h_surface = Nu * k_gas / d_o.

We also compute an effective overall coefficient h_effective referenced to
the total external area:

  UA_effective = h_surface * (A_tube_total + η_f * A_fin_total)
  h_effective  = UA_effective / A_total

where η_f is the fin efficiency. If the user does not provide a value,
we estimate η_f using a Schmidt-type approximation for circular fins
based on the local h_surface and a fin thermal conductivity k_fin.
Finally, we report the classical fixed-pumping-power performance index:

  PEC = j / f**(1/3).
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import warnings
from typing import Optional


@dataclass
class GasProperties:
    """
    Gas-side thermophysical properties.

    All properties are assumed constant and evaluated at a representative
    film temperature. Default values correspond to air at 0°C (273.15 K).
    """

    rho: float = 1.293  # [kg/m^3]
    cp: float = 1005.9  # [J/kg-K]
    mu: float = 1.7258e-5  # [Pa·s] dynamic viscosity
    k: float = 0.024009  # [W/m-K] thermal conductivity

    @property
    def Pr(self) -> float:
        """Prandtl number Pr = (mu * cp) / k."""
        return self.mu * self.cp / self.k

    @property
    def nu(self) -> float:
        """Kinematic viscosity ν = mu / rho [m^2/s]."""
        return self.mu / self.rho

    @property
    def alpha(self) -> float:
        """Thermal diffusivity α = k / (rho * cp) [m^2/s]."""
        return self.k / (self.rho * self.cp)


@dataclass
class BankGeometry:
    """
    기하 변수 정의 (tube bank 외형 및 튜브/휜 형상).

    All lengths in [m].

    - W_bank: [m] Tube bank width (frontal width, normal to flow).
    - H_bank: [m] Tube bank height (tube axis direction, fin stack direction).
    - L_bank: [m] Tube bank length in flow direction (depth).

    - d_o   : [m] Tube outer diameter.
    - S_T   : [m] Transverse tube pitch (S_t).
    - S_L   : [m] Longitudinal tube pitch (S_l).

    - S_fin : [m] Clear fin spacing s between adjacent fins (Briggs s).
              Fin pitch (center-to-center) = t_fin + S_fin.
    - t_fin : [m] Fin thickness (Briggs t).
    - h_fin : [m] Radial fin height (Briggs h_f).
    """

    W_bank: float
    H_bank: float
    L_bank: float
    d_o: float
    S_T: float
    S_L: float
    S_fin: float
    t_fin: float
    h_fin: float


@dataclass
class BankPerformance:
    """
    결과 값을 담는 데이터 클래스.
    """

    # Flow & dimensionless groups
    Re_d: float
    Nu: float
    j: float
    f: float
    j_over_f_cubert: float

    # Velocities, flow rate, pressure drop, pumping power
    V_fr: float
    V_max: float
    m_dot: float
    delta_p: float
    P_pump: float

    # Geometry counts
    N_tubes: int
    N_rows: int
    N_cols: int
    N_fins_per_tube: int
    N_fins_total: int

    # Areas
    A_fr: float
    A_min: float
    sigma: float
    A_fin_total: float
    A_tube_total: float
    A_total: float

    # Heat-transfer coefficients
    h_surface: float
    h_effective: float
    UA_effective: float


class StaggeredFinnedTubeBank:
    """
    Staggered circular finned-tube bank model using Briggs & Young (1963)
    heat transfer correlation and ESDU 86022 (1986) pressure drop correlation.

    Typical usage
    -------------
    >>> gas = GasProperties()  # air at 0°C
    >>> geom = BankGeometry(
    ...     W_bank=1.0,
    ...     H_bank=0.5,
    ...     L_bank=1.0,
    ...     d_o=0.025,
    ...     S_T=0.04,
    ...     S_L=0.035,
    ...     S_fin=0.002,
    ...     t_fin=0.0005,
    ...     h_fin=0.01,
    ... )
    >>> bank = StaggeredFinnedTubeBank(
    ...     geometry=geom,
    ...     gas=gas,
    ...     P_pump_target=3000.0,
    ...     fin_efficiency=None,
    ...     k_fin=200.0,
    ... )
    >>> perf = bank.evaluate(T_s=40.0, T_g_in=10.0)
    >>> perf.j_over_f_cubert  # screening metric under fixed pumping power
    0.5  # doctest: +SKIP
    """

    def __init__(
        self,
        geometry: BankGeometry,
        gas: Optional[GasProperties] = None,
        P_pump_target: float = 3000.0,
        fin_efficiency: Optional[float] = None,
        k_fin: float = 200.0,
    ) -> None:
        """
        Parameters
        ----------
        geometry:
            BankGeometry instance with external volume and finned-tube layout.
        gas:
            GasProperties instance. If None, defaults to air at 0°C.
        P_pump_target:
            Target pumping power on gas side [W].
        fin_efficiency:
            Overall fin efficiency η_f (0–1).
            - If a float is provided, it is used directly.
            - If None, η_f is estimated using a Schmidt-type approximation
              for circular fins based on the local h_surface and k_fin.
        k_fin:
            [W/m-K] Thermal conductivity of fin material (e.g. 200 W/m-K
            for aluminum).
        """
        self.geometry = geometry
        self.gas = gas if gas is not None else GasProperties()
        self.P_pump_target = float(P_pump_target)
        self._fin_efficiency_input = fin_efficiency
        self.k_fin = float(k_fin)

        if self.geometry.S_T <= self.geometry.d_o:
            raise ValueError(
                "Transverse pitch S_T must be larger than tube outer diameter d_o "
                "to have a positive minimum flow area."
            )

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    def _compute_tube_and_fin_counts(self) -> tuple[int, int, int, int, int]:
        """
        Compute tube counts and fin counts based on bank external volume
        and pitch values.

        Assumptions
        -----------
        - Tubes are arranged in a staggered array in the frontal plane.
        - N_cols (number of tubes across the width W_bank) is approximated as:
              N_cols = floor(W_bank / S_T)
        - N_rows (number of tube rows in the flow direction) is:
              N_rows = floor(L_bank / S_L)
        - Each tube extends over the full bank height H_bank, and fins are
          stacked uniformly along this height.
        - Fin pitch (center-to-center) = t_fin + S_fin, so
              N_fins_per_tube = floor(H_bank / (t_fin + S_fin)).
        """
        g = self.geometry

        N_cols = max(1, int(math.floor(g.W_bank / g.S_T)))
        N_rows = max(1, int(math.floor(g.L_bank / g.S_L)))
        N_tubes = N_cols * N_rows

        fin_pitch = g.t_fin + g.S_fin
        if fin_pitch <= 0.0:
            raise ValueError("Fin pitch (t_fin + S_fin) must be positive.")

        N_fins_per_tube = max(1, int(math.floor(g.H_bank / fin_pitch)))
        N_fins_total = N_fins_per_tube * N_tubes

        return N_tubes, N_rows, N_cols, N_fins_per_tube, N_fins_total

    def _compute_flow_areas(self) -> tuple[float, float, float]:
        """
        Compute frontal area A_fr, minimum flow area A_min, and area ratio σ.

        Definitions and approximations
        ------------------------------
        - Frontal area:
              A_fr = W_bank * H_bank

        - Minimum flow area A_min:
          For a staggered tube bank, the minimum free-flow area can occur
          either between tubes in the transverse direction or along a
          diagonal path between staggered tubes. To reflect this, we define
          an effective gap width g_min as the smaller of:

              g_T = S_T - d_o                  (transverse gap)
              g_D = S_d - d_o,  S_d = sqrt( (S_T/2)^2 + S_L^2 )

          and set

              σ = A_min / A_fr ≈ g_min / S_T
              A_min = σ * A_fr

          This keeps the Reynolds number definition consistent with
          Briggs & Young practice while accounting for the possibility
          that the diagonal clearance controls the maximum velocity in
          a staggered array.

          Note: Effects of thin circular fins on the minimum area are
          handled implicitly via the heat-transfer correlation; we use a
          bare-tube-style A_min for Re_d.
        """
        g = self.geometry

        A_fr = g.W_bank * g.H_bank

        S_t = g.S_T
        S_l = g.S_L
        d_o = g.d_o

        g_T = S_t - d_o
        if g_T <= 0.0:
            raise ValueError(
                "Transverse gap (S_T - d_o) must be positive to obtain "
                "a finite minimum flow area."
            )

        S_d = math.sqrt((S_t / 2.0) ** 2 + S_l**2)
        g_D = S_d - d_o

        if g_D <= 0.0:
            warnings.warn(
                "Diagonal gap g_D = S_d - d_o is non-positive; using "
                "transverse gap only for A_min.",
                RuntimeWarning,
            )
            g_min = g_T
        else:
            g_min = min(g_T, g_D)

        sigma = g_min / S_t
        if sigma <= 0.0:
            raise ValueError(
                "Computed area ratio sigma <= 0; check S_T, S_L, and d_o values."
            )

        A_min = A_fr * sigma
        return A_fr, A_min, sigma

    # ------------------------------------------------------------------
    # Correlation implementations
    # ------------------------------------------------------------------

    def _briggs_young_heat_transfer(self, Re_d: float) -> tuple[float, float]:
        """
        Briggs & Young (1963) Nusselt number and Colburn j-factor.

        Nu = 0.134 Re_d^0.681 Pr^(1/3) (s/h_f)^0.2 (s/t)^0.1134

        Returns
        -------
        Nu, j
        """
        g = self.geometry
        gp = self.gas

        s = g.S_fin
        h_f = g.h_fin
        t = g.t_fin

        if h_f <= 0.0:
            raise ValueError("Fin height h_fin (h_f) must be positive.")
        if s <= 0.0:
            raise ValueError("Clear fin spacing S_fin (s) must be positive.")
        if t <= 0.0:
            raise ValueError("Fin thickness t_fin (t) must be positive.")

        Nu = (
            0.134
            * (Re_d ** 0.681)
            * (gp.Pr ** (1.0 / 3.0))
            * (s / h_f) ** 0.2
            * (s / t) ** 0.1134
        )
        j = Nu / (Re_d * (gp.Pr ** (1.0 / 3.0)))

        # Validity check: Re_d ~ 1.1e3 – 1.8e4 for Briggs & Young
        if not (1.1e3 <= Re_d <= 1.8e4):
            warnings.warn(
                f"Re_d={Re_d:.1f} is outside the nominal range "
                "of the Briggs & Young correlation (1.1e3–1.8e4).",
                RuntimeWarning,
            )

        return Nu, j

    def _esdu_86022_friction(
        self,
        Re_d: float,
        A_increase: float,
        sigma: float,
    ) -> tuple[float, float]:
        """
        ESDU 86022 (1986) pressure drop correlation for high-fin staggered
        tube banks.

        K_f = 4.567 * Re_d^(-0.242) * A_increase^0.504
              * (S_T/d_o)^(-0.376) * (S_L/d_o)^(-0.546)

        K_acc = 1 + σ²

        ΔP = (K_acc + N_rows * K_f) * (1/2) * ρ * V_max²

        Returns
        -------
        K_f : float
            Friction factor per tube row.
        K_acc : float
            Acceleration loss coefficient.
        """
        g = self.geometry

        S_T = g.S_T
        S_L = g.S_L
        d_o = g.d_o

        # ESDU 86022 friction factor
        K_f = (
            4.567
            * (Re_d ** -0.242)
            * (A_increase ** 0.504)
            * (S_T / d_o) ** -0.376
            * (S_L / d_o) ** -0.546
        )

        # Acceleration loss
        K_acc = 1.0 + sigma * sigma

        # Validity check for ESDU 86022
        # Re 범위: 10² ~ 10⁵ (Gemini 문헌 조사 결과)
        # 튜브 피치: 1.1 ≤ S_T/d_o ≤ 4.0, 1.1 ≤ S_L/d_o ≤ 3.0
        if not (1.0e2 <= Re_d <= 1.0e5):
            warnings.warn(
                f"Re_d={Re_d:.1f} is outside the ESDU 86022 validated range "
                "(10² ~ 10⁵).",
                RuntimeWarning,
            )

        ST_ratio = S_T / d_o
        SL_ratio = S_L / d_o
        if not (1.1 <= ST_ratio <= 4.0):
            warnings.warn(
                f"S_T/d_o={ST_ratio:.2f} is outside the ESDU 86022 validated range "
                "(1.1 ~ 4.0).",
                RuntimeWarning,
            )
        if not (1.1 <= SL_ratio <= 3.0):
            warnings.warn(
                f"S_L/d_o={SL_ratio:.2f} is outside the ESDU 86022 validated range "
                "(1.1 ~ 3.0).",
                RuntimeWarning,
            )

        return K_f, K_acc

    def _compute_A_increase(
        self,
        N_tubes: int,
        N_fins_per_tube: int,
    ) -> float:
        """
        Compute A_increase = A_total / A_bare_tube for ESDU 86022.

        A_total: 전체 핀+튜브 표면적
        A_bare_tube: 핀이 없다고 가정했을 때의 베어튜브 외부 표면적
        """
        g = self.geometry

        # 핀 면적 (단일 핀)
        r_i = 0.5 * g.d_o
        r_o = r_i + g.h_fin
        A_fin_single = (
            2.0 * math.pi * (r_o**2 - r_i**2)
            + 2.0 * math.pi * r_o * g.t_fin
        )

        N_fins_total = N_tubes * N_fins_per_tube
        A_fin_total = A_fin_single * N_fins_total

        # 튜브 표면적 (핀 사이 노출 부분)
        A_tube_per_tube = math.pi * g.d_o * g.H_bank
        A_tube_total = A_tube_per_tube * N_tubes

        A_total = A_fin_total + A_tube_total

        # 베어튜브 표면적 (핀이 없다고 가정)
        A_bare_tube = A_tube_per_tube * N_tubes

        if A_bare_tube <= 0.0:
            return 1.0

        return A_total / A_bare_tube

    # ------------------------------------------------------------------
    # Pumping power constraint and root finding
    # ------------------------------------------------------------------

    def _pumping_power_from_Vfr(
        self,
        V_fr: float,
        A_fr: float,
        sigma: float,
        N_rows: int,
        A_increase: float,
    ) -> tuple[float, float, float, float]:
        """
        Compute pumping power and related quantities for a given frontal
        velocity V_fr using ESDU 86022 pressure drop correlation.

        Returns
        -------
        P_pump, Re_d, K_f, delta_p

        Notes
        -----
        K_f는 ESDU 86022의 마찰계수 (단일 열 기준)로, Robinson & Briggs의
        Fanning friction factor f와는 다른 정의를 가진다.
        """
        gp = self.gas
        g = self.geometry

        V_max = V_fr / sigma
        Re_d = gp.rho * V_max * g.d_o / gp.mu

        # ESDU 86022 pressure drop
        K_f, K_acc = self._esdu_86022_friction(
            Re_d=Re_d,
            A_increase=A_increase,
            sigma=sigma,
        )

        # Pressure drop: ΔP = (K_acc + N_rows * K_f) * (1/2) * ρ * V_max²
        delta_p = (K_acc + N_rows * K_f) * (gp.rho * V_max**2 / 2.0)

        # Mass flow rate and pumping power
        m_dot = gp.rho * V_fr * A_fr
        P_pump = (m_dot / gp.rho) * delta_p

        return P_pump, Re_d, K_f, delta_p

    def _solve_velocity_for_pumping_power(
        self,
        A_fr: float,
        sigma: float,
        N_rows: int,
        A_increase: float,
        V_fr_min: float = 0.01,
        V_fr_max_initial: float = 50.0,
        tol_rel: float = 1e-6,
        max_iter: int = 60,
    ) -> tuple[float, float, float, float]:
        """
        Find the frontal velocity V_fr that satisfies

            P_pump(V_fr) = P_pump_target

        using a bisection method with adaptive upper bound expansion.
        """
        target = self.P_pump_target

        # Handle trivial zero target
        if target <= 0.0:
            return 0.0, 0.0, 0.0, 0.0

        # Ensure lower bound
        V_lo = max(V_fr_min, 1e-6)
        P_lo, Re_lo, f_lo, dp_lo = self._pumping_power_from_Vfr(
            V_lo, A_fr, sigma, N_rows, A_increase
        )

        # If even at very small V_fr power already exceeds target, warn
        if P_lo > target:
            warnings.warn(
                "Even the minimum frontal velocity considered yields "
                "pumping power above the target. Using V_fr = V_fr_min.",
                RuntimeWarning,
            )
            return V_lo, Re_lo, f_lo, dp_lo

        # Find upper bound by expansion
        V_hi = V_fr_max_initial
        P_hi, Re_hi, f_hi, dp_hi = self._pumping_power_from_Vfr(
            V_hi, A_fr, sigma, N_rows, A_increase
        )
        expansion_count = 0
        while P_hi < target and expansion_count < 10:
            V_hi *= 2.0
            P_hi, Re_hi, f_hi, dp_hi = self._pumping_power_from_Vfr(
                V_hi, A_fr, sigma, N_rows, A_increase
            )
            expansion_count += 1

        if P_hi < target:
            warnings.warn(
                "Unable to bracket the target pumping power even after "
                "expanding V_fr; using V_fr = V_hi.",
                RuntimeWarning,
            )
            return V_hi, Re_hi, f_hi, dp_hi

        # Bisection
        for _ in range(max_iter):
            V_mid = 0.5 * (V_lo + V_hi)
            P_mid, Re_mid, f_mid, dp_mid = self._pumping_power_from_Vfr(
                V_mid, A_fr, sigma, N_rows, A_increase
            )

            if abs(P_mid - target) <= tol_rel * target:
                return V_mid, Re_mid, f_mid, dp_mid

            if P_mid < target:
                V_lo = V_mid
                P_lo = P_mid
            else:
                V_hi = V_mid
                P_hi = P_mid

        # If not converged within max_iter, return last midpoint
        warnings.warn(
            "Bisection did not fully converge within max_iter; "
            "returning last midpoint.",
            RuntimeWarning,
        )
        return V_mid, Re_mid, f_mid, dp_mid

    # ------------------------------------------------------------------
    # Area and heat-transfer coefficient calculations
    # ------------------------------------------------------------------

    def _compute_areas(
        self,
        N_tubes: int,
        N_fins_per_tube: int,
    ) -> tuple[float, float, float]:
        """
        Compute fin, tube, and total external areas.

        Fin area (one fin):
          - Outer radius r_o = d_o/2 + h_fin
          - Inner radius r_i = d_o/2
          - Disk faces (both sides): 2 * π * (r_o^2 - r_i^2)
          - Outer perimeter edge:   2 * π * r_o * t_fin

        A_fin_single = 2π(r_o^2 - r_i^2) + 2π r_o t_fin
        """
        g = self.geometry

        r_i = 0.5 * g.d_o
        r_o = r_i + g.h_fin

        A_fin_single = (
            2.0 * math.pi * (r_o**2 - r_i**2)
            + 2.0 * math.pi * r_o * g.t_fin
        )

        N_fins_total = N_tubes * N_fins_per_tube
        A_fin_total = A_fin_single * N_fins_total

        # Approximate total tube external area (ignoring fin shading)
        A_tube_per_tube = math.pi * g.d_o * g.H_bank
        A_tube_total = A_tube_per_tube * N_tubes

        A_total = A_fin_total + A_tube_total
        return A_fin_total, A_tube_total, A_total

    def _compute_fin_efficiency(self, h_surface: float) -> float:
        """
        Schmidt-type approximate efficiency for circular fins.

        We approximate the circular (annular) fin as an equivalent straight
        fin of length L = h_fin, thickness t_fin, and uniform convection
        coefficient h_surface, leading to

            m = sqrt(2 h / (k_fin t_fin))
            η_f ≈ tanh(m L) / (m L)

        This is a common engineering approximation for high-fin circular
        fins when detailed Bessel-function formulas are not required.
        """
        g = self.geometry

        if g.t_fin <= 0.0 or g.h_fin <= 0.0:
            return 1.0

        m = math.sqrt(max(0.0, 2.0 * h_surface / (self.k_fin * g.t_fin)))
        L = g.h_fin

        x = m * L
        if x < 1e-8:
            return 1.0

        eta_f = math.tanh(x) / x
        return max(0.0, min(1.0, eta_f))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def evaluate(
        self,
        T_s: float,
        T_g_in: float,
    ) -> BankPerformance:
        """
        Evaluate bank performance under the fixed pumping power constraint.

        Parameters
        ----------
        T_s:
            [°C] Tube surface temperature. (현재 구현에서는 물성 보정에
            직접 사용하지 않고, 향후 film temperature 기반 물성 업데이트를
            위해 인터페이스만 준비한다.)
        T_g_in:
            [°C] Inlet gas temperature. (위와 동일하게 인터페이스용 인자.)

        Returns
        -------
        BankPerformance
            j, f, j/f^(1/3), Nu, Re_d, ΔP, m_dot, V_fr, V_max,
            fin/tube/total areas, h_surface, h_effective, UA_effective 등을 포함.
        """
        # Current implementation keeps gas properties fixed (e.g. at 0°C).
        # If desired, the user can instantiate GasProperties with film-
        # temperature-based values and pass it to this class.
        _ = (T_s, T_g_in)  # unused but kept for future extension

        # Geometry and counts
        N_tubes, N_rows, N_cols, N_fins_per_tube, N_fins_total = (
            self._compute_tube_and_fin_counts()
        )
        A_fr, A_min, sigma = self._compute_flow_areas()

        # Compute A_increase for ESDU 86022
        A_increase = self._compute_A_increase(
            N_tubes=N_tubes,
            N_fins_per_tube=N_fins_per_tube,
        )

        # Solve for V_fr under pumping power constraint (using ESDU 86022)
        V_fr, Re_d, K_f, delta_p = self._solve_velocity_for_pumping_power(
            A_fr=A_fr,
            sigma=sigma,
            N_rows=N_rows,
            A_increase=A_increase,
        )

        # Dimensionless groups and coefficients at the converged Re_d
        Nu, j = self._briggs_young_heat_transfer(Re_d)
        gp = self.gas
        g = self.geometry

        V_max = V_fr / sigma if sigma > 0.0 else 0.0
        m_dot = gp.rho * V_fr * A_fr
        P_pump = (m_dot / gp.rho) * delta_p

        # NOTE: K_f는 ESDU 86022의 마찰계수 (단일 열 기준).
        # j/f^(1/3)의 f에는 K_f를 사용한다.
        f = K_f
        j_over_f_cubert = j / (f ** (1.0 / 3.0)) if f > 0.0 else float("nan")

        # Areas and effective coefficients
        A_fin_total, A_tube_total, A_total = self._compute_areas(
            N_tubes=N_tubes,
            N_fins_per_tube=N_fins_per_tube,
        )

        h_surface = Nu * gp.k / g.d_o

        if self._fin_efficiency_input is not None:
            eta_f = self._fin_efficiency_input
        else:
            eta_f = self._compute_fin_efficiency(h_surface)

        UA_effective = h_surface * (
            A_tube_total + eta_f * A_fin_total
        )
        h_effective = UA_effective / A_total if A_total > 0.0 else 0.0

        return BankPerformance(
            Re_d=Re_d,
            Nu=Nu,
            j=j,
            f=f,
            j_over_f_cubert=j_over_f_cubert,
            V_fr=V_fr,
            V_max=V_max,
            m_dot=m_dot,
            delta_p=delta_p,
            P_pump=P_pump,
            N_tubes=N_tubes,
            N_rows=N_rows,
            N_cols=N_cols,
            N_fins_per_tube=N_fins_per_tube,
            N_fins_total=N_fins_total,
            A_fr=A_fr,
            A_min=A_min,
            sigma=sigma,
            A_fin_total=A_fin_total,
            A_tube_total=A_tube_total,
            A_total=A_total,
            h_surface=h_surface,
            h_effective=h_effective,
            UA_effective=UA_effective,
        )


if __name__ == "__main__":
    # 간단한 예시 실행용 스니펫.
    gas = GasProperties()
    geom = BankGeometry(
        W_bank=0.42,
        H_bank=0.42,
        L_bank=0.165,
        d_o=0.01588,
        S_T=0.04,
        S_L=0.035,
        S_fin=0.004,   # ESDU 86022 범위 내 (s/d_o ≈ 0.25)
        t_fin=0.00025,
        h_fin=0.007,
    )
    bank = StaggeredFinnedTubeBank(
        geometry=geom,
        gas=gas,
        P_pump_target=1000.0,
        fin_efficiency=None,  # use Schmidt approximation
        k_fin=200.0,  # [W/m-K], e.g. aluminum
    )
    perf = bank.evaluate(T_s=40.0, T_g_in=10.0)

    print("=== ESDU 86022 Pressure Drop Test ===")
    print("Re_d          :", f"{perf.Re_d:.1f}")
    print("j             :", f"{perf.j:.5f}")
    print("K_f (ESDU)    :", f"{perf.f:.5f}")
    print("j/f^(1/3)     :", f"{perf.j_over_f_cubert:.5f}")
    print("A_fin_total   :", f"{perf.A_fin_total:.4f}", "[m^2]")
    print("A_tube_total  :", f"{perf.A_tube_total:.4f}", "[m^2]")
    print("A_total       :", f"{perf.A_total:.4f}", "[m^2]")
    print("h_surface     :", f"{perf.h_surface:.2f}", "[W/m^2-K]")
    print("h_effective   :", f"{perf.h_effective:.2f}", "[W/m^2-K]")
    print("UA_effective  :", f"{perf.UA_effective:.2f}", "[W/K]")
