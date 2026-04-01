"""Generate fig_result_shared_fp: Correct case (FP trigger) vs gravity-flipped error."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata
from pathlib import Path
import os

for font in ['Arial', 'Liberation Sans', 'DejaVu Sans']:
    try:
        matplotlib.font_manager.findfont(font, fallback_to_default=False)
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = [font]
        FONT_USED = font
        break
    except:
        continue
else:
    FONT_USED = 'DejaVu Sans'

plt.rcParams['font.size'] = 7

BASE = Path(__file__).resolve().parent.parent
CASES = BASE / 'benchmark' / 'cases'
OUT_DIR = BASE / 'results'


def load_velocity(case_dir, vtm_name):
    vtm_path = case_dir / 'VTK' / vtm_name
    mesh = pv.read(vtm_path)
    block = mesh[0]
    points = block.points
    U = block.point_data['U']
    mag = np.linalg.norm(U, axis=1)
    return points, mag


def interpolate_to_grid(points, vals, nx=800, ny=300):
    x, y = points[:, 0], points[:, 1]
    x0, x1 = float(np.floor(x.min())), float(np.ceil(x.max()))
    y0, y1 = float(np.floor(y.min())), float(np.ceil(y.max()))
    xi = np.arange(nx) * (x1 - x0) / (nx - 1) + x0
    yi = np.arange(ny) * (y1 - y0) / (ny - 1) + y0
    Xi, Yi = np.meshgrid(xi, yi)
    Vi = griddata((x, y), vals, (Xi, Yi), method='cubic')
    return xi, yi, Xi, Yi, Vi


def load_velocity_field(case_dir, vtm_name):
    """Load full velocity vector field for streamlines."""
    vtm_path = case_dir / 'VTK' / vtm_name
    mesh = pv.read(vtm_path)
    block = mesh[0]
    pts = block.points
    U = block.point_data['U']
    return pts, U


pts1, U1 = load_velocity_field(CASES / 'S8_correct_lam', 'S8_correct_lam_10000.vtm')
pts2, U2 = load_velocity_field(CASES / 'S8_E8_gravity_flipped', 'S8_E8_gravity_flipped_10000.vtm')

mag1 = np.linalg.norm(U1, axis=1)
mag2 = np.linalg.norm(U2, axis=1)

xi1, yi1, Xi1, Yi1, Mi1 = interpolate_to_grid(pts1, mag1)
xi2, yi2, Xi2, Yi2, Mi2 = interpolate_to_grid(pts2, mag2)

# Interpolate velocity components for streamlines
_, _, _, _, Ux1 = interpolate_to_grid(pts1, U1[:, 0])
_, _, _, _, Uy1 = interpolate_to_grid(pts1, U1[:, 1])
_, _, _, _, Ux2 = interpolate_to_grid(pts2, U2[:, 0])
_, _, _, _, Uy2 = interpolate_to_grid(pts2, U2[:, 1])

vmin = min(mag1.min(), mag2.min())
vmax = max(mag1.max(), mag2.max())

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.5, 3.2),
    gridspec_kw={'wspace': 0.20, 'left': 0.08, 'right': 0.90, 'top': 0.97, 'bottom': 0.12})

# (a) Correct case — triggers false positives
im = ax_a.pcolormesh(Xi1, Yi1, Mi1, cmap='viridis', vmin=vmin, vmax=vmax,
                     shading='auto', rasterized=True)
# Replace NaN with 0 for streamplot
Ux1_s = np.nan_to_num(Ux1, nan=0.0)
Uy1_s = np.nan_to_num(Uy1, nan=0.0)
ax_a.streamplot(xi1, yi1, Ux1_s, Uy1_s, color='white', linewidth=0.3,
                density=1.5, arrowsize=0.4, arrowstyle='->')
ax_a.set_xlabel('x [m]', fontsize=7)
ax_a.set_ylabel('y [m]', fontsize=7)
ax_a.set_aspect('auto')
ax_a.tick_params(labelsize=6)
ax_a.text(0.03, 0.92, '(a)', transform=ax_a.transAxes, fontsize=7,
          fontweight='semibold', va='top')

# (b) Gravity-flipped error
im_b = ax_b.pcolormesh(Xi2, Yi2, Mi2, cmap='viridis', vmin=vmin, vmax=vmax,
                       shading='auto', rasterized=True)
Ux2_s = np.nan_to_num(Ux2, nan=0.0)
Uy2_s = np.nan_to_num(Uy2, nan=0.0)
ax_b.streamplot(xi2, yi2, Ux2_s, Uy2_s, color='white', linewidth=0.3,
                density=1.5, arrowsize=0.4, arrowstyle='->')
ax_b.set_xlabel('x [m]', fontsize=7)
ax_b.set_ylabel('')
ax_b.set_aspect('auto')
ax_b.tick_params(labelsize=6)
ax_b.text(0.03, 0.92, '(b)', transform=ax_b.transAxes, fontsize=7,
          fontweight='semibold', va='top')

# Shared colorbar on right side
cb = fig.colorbar(im, ax=[ax_a, ax_b], shrink=0.85, pad=0.015, aspect=20)
cb.set_label('|U| [m/s]', fontsize=7)
cb.ax.tick_params(labelsize=6)

out_pdf = OUT_DIR / 'fig_result_shared_fp.pdf'
out_png = OUT_DIR / 'fig_result_shared_fp.png'
fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
plt.close()

print(f"Font used: {FONT_USED}")
print(f"|U| range: {vmin:.4f} - {vmax:.4f} m/s")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
