"""Generate fig_result_accuracy_and_heatmap: Merged accuracy bar + error heatmap."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

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
OUT_DIR = BASE / 'results'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 2.6),
                                gridspec_kw={'width_ratios': [0.85, 1.15], 'wspace': 0.40})

# ── (a) Accuracy bar chart ──
evaluators = ['Claude', 'Gemini', 'Expert']
setup_acc = [99.6, 87.5, 73.8]
imgonly_acc = [83.3, 79.2, 66.7]

x = np.arange(len(evaluators))
w = 0.30

bars1 = ax1.bar(x - w/2, setup_acc, w, label='Setup-conditioned',
                color='#4A7FC1', edgecolor='white', linewidth=0.5)
bars2 = ax1.bar(x + w/2, imgonly_acc, w, label='Image-only',
                color='#B0C8E8', edgecolor='white', linewidth=0.5)

for bar in bars1:
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
             f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=6)
for bar in bars2:
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
             f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=6)

ax1.set_ylabel('Accuracy (%)', fontsize=7)
ax1.set_ylim(0, 108)
ax1.set_xticks(x)
ax1.set_xticklabels(evaluators, fontsize=7)
ax1.tick_params(labelsize=6)
ax1.legend(fontsize=5.5, loc='upper right', frameon=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.text(-0.18, 1.03, '(a)', transform=ax1.transAxes, fontsize=7, fontweight='semibold')

# ── (b) Error-type recall heatmap ──
error_labels = ['E1', 'E2', 'E3', 'E4', 'E5', 'E8']

recall_data = np.array([
    [100, 100, 100, 100, 100, 100],   # Claude
    [100,  50,  67, 100,  50,  80],   # Gemini
    [ 73,  20,  29,  50,  44, 100],   # Expert
])

eval_labels = ['Claude', 'Gemini', 'Expert']

im = ax2.imshow(recall_data, cmap='viridis', vmin=0, vmax=100, aspect='auto')

ax2.set_xticks(range(len(error_labels)))
ax2.set_xticklabels(error_labels, fontsize=7)
ax2.set_yticks(range(len(eval_labels)))
ax2.set_yticklabels(eval_labels, fontsize=7)

for i in range(len(eval_labels)):
    for j in range(len(error_labels)):
        val = recall_data[i, j]
        # viridis: dark purple (0) → teal (50) → yellow (100)
        # white text readable on dark purple/teal (0-55), black on green/yellow (56+)
        color = 'white' if val <= 55 else 'black'
        ax2.text(j, i, f'{val}', ha='center', va='center',
                 fontsize=7, color=color)

ax2.set_xlabel('Error type', fontsize=7)
ax2.text(-0.15, 1.03, '(b)', transform=ax2.transAxes, fontsize=7, fontweight='semibold')

cbar = fig.colorbar(im, ax=ax2, shrink=0.85, pad=0.03)
cbar.set_label('Recall (%)', fontsize=7)
cbar.ax.tick_params(labelsize=6)

plt.tight_layout()

out_pdf = OUT_DIR / 'fig_result_accuracy_and_heatmap.pdf'
out_png = OUT_DIR / 'fig_result_accuracy_and_heatmap.png'
fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
plt.close()

print(f"Font used: {FONT_USED}")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
