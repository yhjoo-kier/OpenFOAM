"""Fig 7: (a) Per-trial accuracy strip plot, (b) contamination correction."""
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

models = ['Claude', 'GPT-5.4', 'Gemini']
colors = ['#D4770C', '#10A37F', '#4285F4']

setup_trials = {
    'Claude': [90.0, 90.0, 86.7],
    'GPT-5.4': [80.0, 80.0, 80.0],
    'Gemini': [56.7, 53.3, 63.3],
}
io_trials = {
    'Claude': [33.3, 33.3, 33.3],
    'GPT-5.4': [43.3, 40.0, 46.7],
    'Gemini': [63.3, 56.7, 40.0],
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.6),
    gridspec_kw={'width_ratios': [1.2, 0.8], 'wspace': 0.35})

# ── (a) Strip plot with mean + individual trials ──
y_base = {'Claude': 0, 'GPT-5.4': 1, 'Gemini': 2}
trial_jitter = [-0.06, 0, 0.06]

for mi, (model, c) in enumerate(zip(models, colors)):
    y = y_base[model]
    s_vals = setup_trials[model]
    io_vals = io_trials[model]
    s_mean = np.mean(s_vals)
    io_mean = np.mean(io_vals)

    # Setup trials (filled circles, upper row)
    for t, val in enumerate(s_vals):
        ax1.scatter(val, y - 0.18 + trial_jitter[t], color=c, s=20, zorder=3,
                    marker='o', edgecolors='none', alpha=0.7)
    # Setup mean (larger filled diamond)
    ax1.scatter(s_mean, y - 0.18, color=c, s=50, zorder=4, marker='D',
                edgecolors='white', linewidths=0.5)

    # IO trials (open circles, lower row)
    for t, val in enumerate(io_vals):
        ax1.scatter(val, y + 0.18 + trial_jitter[t], color=c, s=20, zorder=3,
                    marker='o', facecolors='none', edgecolors=c, linewidths=0.8, alpha=0.7)
    # IO mean (larger open diamond)
    ax1.scatter(io_mean, y + 0.18, color=c, s=50, zorder=4, marker='D',
                facecolors='none', edgecolors=c, linewidths=1.0)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=4, label='Trial (setup-cond.)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor='gray',
           markersize=4, markeredgewidth=0.8, label='Trial (image-only)'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='gray', markersize=5, label='Mean'),
]
ax1.legend(handles=legend_elements, fontsize=5, loc='lower left', frameon=False)

ax1.set_yticks([0, 1, 2])
ax1.set_yticklabels(models, fontsize=6)
ax1.set_xlabel('Accuracy (%)', fontsize=7)
ax1.set_xlim(25, 95)
ax1.tick_params(labelsize=6)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.grid(True, alpha=0.2, linewidth=0.5)
ax1.set_axisbelow(True)
ax1.text(-0.15, 1.03, '(a)', transform=ax1.transAxes, fontsize=7, fontweight='semibold')

# ── (b) Contamination correction ──
labels = ['Subagent\n(contaminated)', 'CLI\n(partial isol.)', 'API\n(full isol.)']
values = [99.6, 83.3, 88.9]
bar_colors = ['#FFCDD2', '#FFF9C4', '#C8E6C9']
edge_colors = ['#C62828', '#F9A825', '#2E7D32']

bars = ax2.barh(range(len(labels)), values, color=bar_colors, edgecolor=edge_colors,
                linewidth=0.8, height=0.5)

for j, (v, ec) in enumerate(zip(values, edge_colors)):
    ax2.text(v + 0.8, j, f'{v:.1f}%', va='center', ha='left', fontsize=6, color=ec)

ax2.set_yticks(range(len(labels)))
ax2.set_yticklabels(labels, fontsize=6)
ax2.set_xlabel('Claude setup-cond. accuracy (%)', fontsize=6)
ax2.set_xlim(0, 108)
ax2.tick_params(labelsize=6)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.text(-0.22, 1.03, '(b)', transform=ax2.transAxes, fontsize=7, fontweight='semibold')

fig.savefig(OUT / 'fig_result_trial_consistency.pdf', format='pdf', bbox_inches='tight')
fig.savefig(OUT / 'fig_result_trial_consistency.png', format='png', dpi=600, bbox_inches='tight')
plt.close()
print("Saved: fig_result_trial_consistency.pdf/.png")
