"""Generate fig_result_bc_swap_illusion: re-rendered from VTK with shared colorbar."""
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

plt.rcParams['font.size'] = 8

BASE = Path(__file__).resolve().parent.parent
CASES = BASE / 'benchmark' / 'cases'
OUT_DIR = BASE / 'results'


def load_temperature(case_dir):
    vtk_dir = case_dir / 'VTK'
    vtm_files = [f for f in os.listdir(vtk_dir) if f.endswith('.vtm') and 'series' not in f]
    mesh = pv.read(vtk_dir / vtm_files[0])
    block = mesh[0]
    points = block.points
    T = block.point_data['T']
    return points, T


def interpolate_to_grid(points, T, nx=200, ny=200):
    x, y = points[:, 0], points[:, 1]
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    Xi, Yi = np.meshgrid(xi, yi)
    Ti = griddata((x, y), T, (Xi, Yi), method='cubic')
    return xi, yi, Ti


pts1, T1 = load_temperature(CASES / 'S6_correct_Ra1e4')
pts2, T2 = load_temperature(CASES / 'S6_E2_bc_swap')

xi1, yi1, Ti1 = interpolate_to_grid(pts1, T1)
xi2, yi2, Ti2 = interpolate_to_grid(pts2, T2)

# Shared colorbar limits
vmin = min(T1.min(), T2.min())
vmax = max(T1.max(), T2.max())

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.5, 3.0),
    gridspec_kw={'wspace': 0.25})

# (a) Correct
im_a = ax_a.pcolormesh(xi1, yi1, Ti1, cmap='coolwarm', vmin=vmin, vmax=vmax,
                        shading='auto', rasterized=True)
ax_a.set_xlabel('x [m]', fontsize=8)
ax_a.set_ylabel('y [m]', fontsize=8)
ax_a.set_aspect('equal')
ax_a.tick_params(labelsize=7)
ax_a.text(0.04, 0.95, '(a)', transform=ax_a.transAxes, fontsize=9,
          fontweight='bold', va='top')

# (b) BC swap
im_b = ax_b.pcolormesh(xi2, yi2, Ti2, cmap='coolwarm', vmin=vmin, vmax=vmax,
                        shading='auto', rasterized=True)
ax_b.set_xlabel('x [m]', fontsize=8)
ax_b.set_ylabel('y [m]', fontsize=8)
ax_b.set_aspect('equal')
ax_b.tick_params(labelsize=7)
ax_b.text(0.04, 0.95, '(b)', transform=ax_b.transAxes, fontsize=9,
          fontweight='bold', va='top')

# Shared colorbar between panels
cb = fig.colorbar(im_a, ax=[ax_a, ax_b], shrink=0.9, pad=0.02, aspect=25)
cb.set_label('T [K]', fontsize=8)
cb.ax.tick_params(labelsize=7)

out_pdf = OUT_DIR / 'fig_result_bc_swap_illusion.pdf'
out_png = OUT_DIR / 'fig_result_bc_swap_illusion.png'
fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
plt.close()

print(f"Font used: {FONT_USED}")
print(f"T range: {vmin:.2f} - {vmax:.2f} K")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
