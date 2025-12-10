"""
Extract CFD results from snappyHexMesh case and compare with correlations.
"""
import os
import math

# Physical properties (air at 300K)
RHO = 1.177      # kg/m³
CP = 1006.0      # J/kg-K
MU = 1.846e-5    # Pa·s
K_AIR = 0.0262   # W/m-K
PR = MU * CP / K_AIR  # Prandtl number

# Geometry parameters
PARAMS = {
    "SL": 0.034,           # 34 mm
    "ST": 0.036,           # 36 mm
    "d_o": 0.016,          # 16 mm
    "d_f": 0.036,          # 36 mm (fin outer diameter)
    "t_fin": 0.0005,       # 0.5 mm
    "P_fin": 0.0045,       # 4.5 mm (fin pitch)
    "h_fin": 0.010,        # 10 mm (fin height = (d_f - d_o)/2)
    "N_fins": 3,           # 3 fins per tube
    "N_tubes": 6,          # 6 tubes total (2 REV x 3 tubes/REV)
    "N_rows": 4,           # 4 tube rows in flow direction
}


def parse_postprocessing_file(filepath):
    """Parse OpenFOAM postProcessing file and return last value."""
    if not os.path.exists(filepath):
        return None
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find last non-comment line
    for line in reversed(lines):
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 2:
                return float(parts[-1])
    return None


def calculate_cfd_results(case_dir):
    """Extract and calculate CFD results from case directory."""
    
    # Read data from time-series file (last line)
    pp_dir = os.path.join(case_dir, "postProcessing/fluid/inletMassFlow/0")
    if not os.path.exists(pp_dir):
        pp_dir = os.path.join(case_dir, "postProcessing/fluid/inletMassFlow")
        if os.path.exists(pp_dir):
            subdirs = [d for d in os.listdir(pp_dir) if d.isdigit()]
            if subdirs:
                pp_dir = os.path.join(pp_dir, max(subdirs, key=int))
    
    if not os.path.exists(pp_dir):
        return None
    
    # Get time from last line of data file
    data_file = os.path.join(pp_dir, "surfaceFieldValue.dat")
    if os.path.exists(data_file):
        with open(data_file, 'r') as f:
            lines = [l.strip() for l in f.readlines() if l.strip() and not l.startswith('#')]
        if lines:
            latest = lines[-1].split()[0]
        else:
            latest = "0"
    else:
        latest = "0"
    print(f"Using time step: {latest}")
    
    # Read postProcessing results (OpenFOAM stores time series in single files under time=0)
    base_pp = os.path.join(case_dir, "postProcessing/fluid")
    
    m_dot = abs(parse_postprocessing_file(
        os.path.join(base_pp, "inletMassFlow/0/surfaceFieldValue.dat")
    ) or 0.0)
    
    T_in = parse_postprocessing_file(
        os.path.join(base_pp, "inletTemperature/0/surfaceFieldValue.dat")
    ) or 300.0
    
    T_out = parse_postprocessing_file(
        os.path.join(base_pp, "outletTemperature/0/surfaceFieldValue.dat")
    ) or 300.0
    
    p_in = parse_postprocessing_file(
        os.path.join(base_pp, "inletPressure/0/surfaceFieldValue.dat")
    ) or 0.0
    
    p_out = parse_postprocessing_file(
        os.path.join(base_pp, "outletPressure/0/surfaceFieldValue.dat")
    ) or 0.0
    
    # Calculate derived quantities
    delta_T = T_out - T_in
    delta_p = p_in - p_out
    Q = m_dot * CP * delta_T  # Heat transfer rate [W]
    
    # Calculate heat transfer coefficient
    # Q = h * A_eff * LMTD
    T_solid = 320.0  # Solid temperature
    LMTD = ((T_solid - T_in) - (T_solid - T_out)) / math.log((T_solid - T_in) / (T_solid - T_out)) if delta_T > 0.01 else (T_solid - (T_in + T_out)/2)
    
    # Effective heat transfer area
    N_tubes = PARAMS["N_tubes"]
    N_fins = PARAMS["N_fins"]
    d_o = PARAMS["d_o"]
    d_f = PARAMS["d_f"]
    t_fin = PARAMS["t_fin"]
    P_fin = PARAMS["P_fin"]
    domain_z = N_fins * P_fin
    
    # Fin area (both sides + edge)
    r_i = d_o / 2
    r_o = d_f / 2
    A_fin_single = 2 * math.pi * (r_o**2 - r_i**2) + 2 * math.pi * r_o * t_fin
    A_fin_total = A_fin_single * N_fins * N_tubes
    
    # Bare tube area (between fins)
    S_fin = P_fin - t_fin  # Clear fin spacing
    A_tube_total = math.pi * d_o * (domain_z - N_fins * t_fin) * N_tubes
    
    A_total = A_fin_total + A_tube_total
    
    # Calculate h from Q and LMTD (assuming fin efficiency = 1 for simplicity)
    h_surface = Q / (A_total * LMTD) if (A_total * LMTD) > 0 else 0.0
    
    # Nusselt number
    Nu = h_surface * d_o / K_AIR
    
    # Reynolds number (based on V_max)
    A_fr = PARAMS["ST"] * domain_z
    sigma = (PARAMS["ST"] - d_o) / PARAMS["ST"]  # Simplified
    V_fr = m_dot / (RHO * A_fr)
    V_max = V_fr / sigma
    Re_d = RHO * V_max * d_o / MU
    
    # Friction factor K_f from pressure drop
    # delta_p = (K_acc + N_rows * K_f) * 0.5 * rho * V_max^2
    K_acc = 1 + sigma**2
    N_rows = PARAMS["N_rows"]
    dyn_press = 0.5 * RHO * V_max**2
    K_total = delta_p / dyn_press if dyn_press > 0 else 0
    K_f = (K_total - K_acc) / N_rows if N_rows > 0 else 0
    
    # Colburn j-factor
    j = Nu / (Re_d * PR**(1/3)) if Re_d > 0 else 0
    
    results = {
        "time_step": latest,
        "m_dot": m_dot,
        "T_in": T_in,
        "T_out": T_out,
        "delta_T": delta_T,
        "p_in": p_in,
        "p_out": p_out,
        "delta_p": delta_p,
        "Q": Q,
        "LMTD": LMTD,
        "A_total": A_total,
        "h_surface": h_surface,
        "Nu": Nu,
        "Re_d": Re_d,
        "V_max": V_max,
        "K_f": K_f,
        "j": j,
    }
    
    return results


def briggs_young_correlation(Re_d, s, h_f, t):
    """Briggs & Young (1963) heat transfer correlation."""
    if Re_d <= 0:
        return 0.0, 0.0
    Nu = 0.134 * Re_d**0.681 * PR**(1/3) * (s/h_f)**0.2 * (s/t)**0.1134
    j = Nu / (Re_d * PR**(1/3))
    return Nu, j


def esdu_86022_friction(Re_d, A_increase, ST_do, SL_do):
    """ESDU 86022 friction factor correlation."""
    K_f = 4.567 * Re_d**(-0.242) * A_increase**0.504 * ST_do**(-0.376) * SL_do**(-0.546)
    return K_f


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    study_dir = os.path.dirname(script_dir)
    case_dir = os.path.join(study_dir, "case_snappy_Re2000")
    
    print("=" * 60)
    print("CFD Results Extraction (snappyHexMesh case)")
    print("=" * 60)
    
    results = calculate_cfd_results(case_dir)
    
    if results is None:
        print("No results found!")
        return
    
    print(f"\n--- CFD Results (Time = {results['time_step']}) ---")
    print(f"Mass flow rate:    {results['m_dot']:.6f} kg/s")
    print(f"Inlet temperature: {results['T_in']:.2f} K")
    print(f"Outlet temperature:{results['T_out']:.2f} K")
    print(f"Temperature rise:  {results['delta_T']:.4f} K")
    print(f"Pressure drop:     {results['delta_p']:.2f} Pa")
    print(f"Heat transfer:     {results['Q']:.2f} W")
    print(f"LMTD:              {results['LMTD']:.2f} K")
    print(f"Total area:        {results['A_total']:.6f} m²")
    print(f"h_surface:         {results['h_surface']:.2f} W/m²-K")
    print(f"Re_d:              {results['Re_d']:.0f}")
    print(f"Nu (CFD):          {results['Nu']:.2f}")
    print(f"K_f (CFD):         {results['K_f']:.4f}")
    print(f"j (CFD):           {results['j']:.6f}")
    
    # Correlation predictions
    print("\n--- Correlation Predictions ---")
    d_o = PARAMS["d_o"]
    s = PARAMS["P_fin"] - PARAMS["t_fin"]  # Clear fin spacing
    h_f = PARAMS["h_fin"]
    t = PARAMS["t_fin"]
    
    Nu_corr, j_corr = briggs_young_correlation(results['Re_d'], s, h_f, t)
    print(f"Nu (Briggs&Young): {Nu_corr:.2f}")
    print(f"j (Briggs&Young):  {j_corr:.6f}")
    
    # A_increase for ESDU
    r_i = d_o / 2
    r_o = PARAMS["d_f"] / 2
    A_fin_single = 2 * math.pi * (r_o**2 - r_i**2) + 2 * math.pi * r_o * t
    A_bare_per_fin = math.pi * d_o * PARAMS["P_fin"]
    A_increase = (A_fin_single + A_bare_per_fin) / A_bare_per_fin
    
    ST_do = PARAMS["ST"] / d_o
    SL_do = PARAMS["SL"] / d_o
    
    K_f_corr = esdu_86022_friction(results['Re_d'], A_increase, ST_do, SL_do)
    print(f"K_f (ESDU 86022):  {K_f_corr:.4f}")
    
    # Comparison
    print("\n--- Comparison ---")
    if Nu_corr > 0:
        Nu_error = (results['Nu'] - Nu_corr) / Nu_corr * 100
        print(f"Nu deviation:      {Nu_error:+.1f}%")
    if K_f_corr > 0 and results['K_f'] > 0:
        Kf_error = (results['K_f'] - K_f_corr) / K_f_corr * 100
        print(f"K_f deviation:     {Kf_error:+.1f}%")


if __name__ == "__main__":
    main()

