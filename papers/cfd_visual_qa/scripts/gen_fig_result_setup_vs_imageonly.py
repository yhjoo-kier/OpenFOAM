"""Generate fig_result_setup_vs_imageonly: contour plots from VTK data.

Row (a): S2_E3_wrong_viscosity — velocity magnitude (channel flow)
Row (b): S1_E3_wrong_viscosity — temperature (heated plate)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata
from pathlib import Path

# ---------------------------------------------------------------------------
# Font setup: Arial > Liberation Sans > DejaVu Sans
# ---------------------------------------------------------------------------
FONT_USED = 'DejaVu Sans'
for font in ['Arial', 'Liberation Sans', 'DejaVu Sans']:
    try:
        matched = fm.findfont(fm.FontProperties(family=font), fallback_to_default=False)
        if matched and font.lower().replace(' ', '') in matched.lower().replace(' ', ''):
            FONT_USED = font
            break
    except (ValueError, Exception):
        continue

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = [FONT_USED, 'DejaVu Sans']
plt.rcParams['font.size'] = 9

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(__file__).resolve().parent.parent
CASE_DIR = BASE / 'benchmark' / 'cases'
OUT_DIR = BASE / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)

VTK_S2 = CASE_DIR / 'S2_E3_wrong_viscosity' / 'VTK' / 'S2_E3_wrong_viscosity_195.vtm'
VTK_S1 = CASE_DIR / 'S1_E3_wrong_viscosity' / 'VTK' / 'S1_E3_wrong_viscosity_10000.vtm'


# ---------------------------------------------------------------------------
# Helper: extract 2D slice from 3D extruded mesh
# ---------------------------------------------------------------------------
def extract_2d(block, field, z_tol=0.01):
    """Return (x, y, values) for the mid-z plane of a thin-extruded mesh."""
    pts = block.points  # (N, 3)
    z_mid = 0.5 * (pts[:, 2].min() + pts[:, 2].max())
    mask = np.abs(pts[:, 2] - z_mid) < z_tol
    # If too few points at exact mid, use all points projected to 2D
    if mask.sum() < 100:
        mask = np.ones(len(pts), dtype=bool)
    x = pts[mask, 0]
    y = pts[mask, 1]
    data = block.point_data[field]
    if data.ndim == 2:
        vals = np.linalg.norm(data[mask], axis=1)
    else:
        vals = data[mask]
    return x, y, vals.astype(float)


def interpolate_grid(x, y, vals, nx=400, ny=None):
    """Interpolate scattered (x,y,vals) onto a uniform grid."""
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    if ny is None:
        aspect = y_range / x_range
        ny = max(20, int(nx * aspect))
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = griddata((x, y), vals, (Xi, Yi), method='linear')
    return Xi, Yi, Zi


# ---------------------------------------------------------------------------
# Load and process data
# ---------------------------------------------------------------------------
mesh_s2 = pv.read(VTK_S2)
block_s2 = mesh_s2[0]
x_s2, y_s2, u_s2 = extract_2d(block_s2, 'U')
Xi_s2, Yi_s2, Zi_s2 = interpolate_grid(x_s2, y_s2, u_s2, nx=600)

mesh_s1 = pv.read(VTK_S1)
block_s1 = mesh_s1[0]
x_s1, y_s1, t_s1 = extract_2d(block_s1, 'T')
Xi_s1, Yi_s1, Zi_s1 = interpolate_grid(x_s1, y_s1, t_s1, nx=300)

# ---------------------------------------------------------------------------
# Figure layout: 2 rows x 1 col, stacked vertically
# Single-column journal width ~3.5 in; two-column span ~7 in
# S2 is wide (10:1), S1 is tall (1:2) → set height ratios accordingly
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(
    2, 1,
    figsize=(7.0, 4.0),
    gridspec_kw=dict(
        hspace=0.40,
        height_ratios=[1.0, 1.5],
        left=0.09, right=0.92, top=0.98, bottom=0.06,
    ),
)

# ---------------------------------------------------------------------------
# Panel (a): S2 — velocity magnitude
# ---------------------------------------------------------------------------
ax = axes[0]
pcm = ax.pcolormesh(
    Xi_s2, Yi_s2, Zi_s2,
    cmap='viridis',
    shading='auto',
    rasterized=True,
    vmin=0.0,
    vmax=np.nanpercentile(Zi_s2, 99),
)
ax.set_aspect('auto')
ax.set_xlabel('x [m]', fontsize=8)
ax.set_ylabel('y [m]', fontsize=8)
ax.tick_params(labelsize=7)
cb = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.015)
cb.set_label('|U| [m/s]', fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.text(0.01, 0.96, '(a)', transform=ax.transAxes,
        fontsize=9, fontweight='bold', va='top', ha='left')

# ---------------------------------------------------------------------------
# Panel (b): S1 — temperature
# ---------------------------------------------------------------------------
ax = axes[1]
pcm = ax.pcolormesh(
    Xi_s1, Yi_s1, Zi_s1,
    cmap='coolwarm',
    shading='auto',
    rasterized=True,
    vmin=np.nanmin(Zi_s1),
    vmax=np.nanmax(Zi_s1),
)
ax.set_aspect('auto')
ax.set_xlabel('x [m]', fontsize=8)
ax.set_ylabel('y [m]', fontsize=8)
ax.tick_params(labelsize=7)
cb = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.025)
cb.set_label('T [K]', fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.text(0.01, 0.98, '(b)', transform=ax.transAxes,
        fontsize=9, fontweight='bold', va='top', ha='left')

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_pdf = OUT_DIR / 'fig_result_setup_vs_imageonly.pdf'
out_png = OUT_DIR / 'fig_result_setup_vs_imageonly.png'
fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
plt.close()

print(f"Font used: {FONT_USED}")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
