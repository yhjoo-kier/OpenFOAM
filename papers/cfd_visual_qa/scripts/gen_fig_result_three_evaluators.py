"""Generate fig_result_three_evaluators: correct vs coarse mesh comparison."""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import griddata
import pyvista as pv

FONT_USED = 'DejaVu Sans'
for font in ['Arial', 'Liberation Sans', 'DejaVu Sans']:
    try:
        matplotlib.font_manager.findfont(font, fallback_to_default=False)
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = [font]
        FONT_USED = font
        break
    except Exception:
        continue

plt.rcParams['font.size'] = 7

BASE = Path(__file__).resolve().parent.parent
CASES = BASE / 'benchmark' / 'cases'
OUT_DIR = BASE / 'results'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_umag(case_dir, vtm_name):
    mesh = pv.read(case_dir / 'VTK' / vtm_name)
    block = mesh[0]
    pts = block.points
    U = block.point_data['U']
    mag = np.linalg.norm(U, axis=1)
    return pts, mag


def interp(pts, mag, nx=200, method='cubic'):
    x, y = pts[:, 0], pts[:, 1]
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), nx)
    Xi, Yi = np.meshgrid(xi, yi)
    Mi = griddata((x, y), mag, (Xi, Yi), method=method)
    return Xi, Yi, Mi


# Load both cases
pts_c, mag_c = load_umag(CASES / 'S5_correct_lam', 'S5_correct_lam_771.vtm')
pts_e, mag_e = load_umag(CASES / 'S5_E5_coarse_mesh', 'S5_E5_coarse_mesh_10000.vtm')

Xi_c, Yi_c, Mi_c = interp(pts_c, mag_c, nx=250, method='cubic')
Xi_e, Yi_e, Mi_e = interp(pts_e, mag_e, nx=250, method='nearest')

# Shared limits
vmin = 0
vmax = max(mag_c.max(), mag_e.max())

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.0, 2.8),
    gridspec_kw={'wspace': 0.25})

# (a) Correct fine mesh
im = ax_a.pcolormesh(Xi_c, Yi_c, Mi_c, cmap='viridis', vmin=vmin, vmax=vmax,
                     shading='auto', rasterized=True)
ax_a.set_xlabel('x [m]', fontsize=7)
ax_a.set_ylabel('y [m]', fontsize=7)
ax_a.set_aspect('equal')
ax_a.tick_params(labelsize=6)
ax_a.text(0.04, 0.94, '(a)', transform=ax_a.transAxes, fontsize=7, fontweight='semibold', va='top')

# (b) Coarse mesh
ax_b.pcolormesh(Xi_e, Yi_e, Mi_e, cmap='viridis', vmin=vmin, vmax=vmax,
                shading='auto', rasterized=True)
ax_b.set_xlabel('x [m]', fontsize=7)
ax_b.set_ylabel('')
ax_b.set_aspect('equal')
ax_b.tick_params(labelsize=6)
ax_b.text(0.04, 0.94, '(b)', transform=ax_b.transAxes, fontsize=7, fontweight='semibold', va='top')

# Shared colorbar
cb = fig.colorbar(im, ax=[ax_a, ax_b], shrink=0.85, pad=0.02, aspect=20)
cb.set_label('|U| [m/s]', fontsize=7)
cb.ax.tick_params(labelsize=6)

out_pdf = OUT_DIR / 'fig_result_three_evaluators.pdf'
out_png = OUT_DIR / 'fig_result_three_evaluators.png'
fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
plt.close()

print(f"Font used: {FONT_USED}")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
