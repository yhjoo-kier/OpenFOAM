"""Fig 6: (a) FP/FN scatter, (b) error-type recall heatmap."""
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
        break
    except:
        continue

plt.rcParams['font.size'] = 7

OUT = Path(__file__).resolve().parent.parent / 'results'

# Mean FP/FN from 3 trials (setup-conditioned)
names = ['Claude', 'GPT-5.4', 'Expert', 'Gemini']
fp_mean = [3.3, 1.7, 2.0, 5.3]
fn_mean = [0.0, 4.3, 19.0, 7.3]
colors = ['#D4770C', '#10A37F', '#4A7FC1', '#4285F4']
markers = ['o', 's', '^', 'D']

# Error-type recall (setup-conditioned, mean across trials)
error_labels = ['E1', 'E2', 'E3', 'E4', 'E5', 'E8']
recall_data = np.array([
    [100, 67, 100, 100, 100, 67],   # Claude (FN=0 overall, but per-type from items)
    [67,  100, 33, 100, 67,  100],  # GPT-5.4
    [73,  20,  29, 50,  44,  100],  # Expert
    [100, 33,  0,  100, 67,  33],   # Gemini
])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.6),
    gridspec_kw={'width_ratios': [0.8, 1.2], 'wspace': 0.40})

# ── (a) FP/FN scatter ──
for i, (n, fp, fn, c, m) in enumerate(zip(names, fp_mean, fn_mean, colors, markers)):
    ax1.scatter(fp, fn, c=c, marker=m, s=60, zorder=3, label=n)
    offset_x = 0.3 if n != 'Expert' else 0.3
    offset_y = 0.5 if n != 'Expert' else -1.0
    ax1.annotate(n, (fp, fn), xytext=(offset_x, offset_y), textcoords='offset points',
                 fontsize=5.5, color=c, va='bottom')

ax1.set_xlabel('Mean false positives', fontsize=7)
ax1.set_ylabel('Mean false negatives', fontsize=7)
ax1.set_xlim(-0.5, 7)
ax1.set_ylim(-1, 22)
ax1.tick_params(labelsize=6)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.axhline(0, color='#cccccc', linewidth=0.5)
ax1.text(-0.18, 1.03, '(a)', transform=ax1.transAxes, fontsize=7, fontweight='semibold')

# ── (b) Recall heatmap ──
im = ax2.imshow(recall_data, cmap='viridis', vmin=0, vmax=100, aspect='auto')
ax2.set_xticks(range(len(error_labels)))
ax2.set_xticklabels(error_labels, fontsize=6)
ax2.set_yticks(range(len(names)))
ax2.set_yticklabels(names, fontsize=6)
ax2.set_xlabel('Error type', fontsize=7)

for i in range(len(names)):
    for j in range(len(error_labels)):
        val = recall_data[i, j]
        color = 'white' if val <= 55 else 'black'
        ax2.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=6, color=color)

ax2.text(-0.18, 1.03, '(b)', transform=ax2.transAxes, fontsize=7, fontweight='semibold')

cb = fig.colorbar(im, ax=ax2, shrink=0.85, pad=0.03)
cb.set_label('Recall (%)', fontsize=7)
cb.ax.tick_params(labelsize=5)

fig.savefig(OUT / 'fig_result_fp_fn_profile.pdf', format='pdf', bbox_inches='tight')
fig.savefig(OUT / 'fig_result_fp_fn_profile.png', format='png', dpi=600, bbox_inches='tight')
plt.close()
print("Saved: fig_result_fp_fn_profile.pdf/.png")
