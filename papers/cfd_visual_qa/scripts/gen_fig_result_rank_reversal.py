"""Fig 5: Rank reversal slopegraph — setup-conditioned vs image-only."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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

evaluators = ['Claude Opus 4.6', 'GPT-5.4', 'Expert', 'Gemini 3.1']
setup =      [88.9,               67.8,      73.8,     57.8]
imageonly =   [33.3,               48.9,      66.7,     53.3]
colors =     ['#D4770C',          '#10A37F',  '#4A7FC1', '#4285F4']
markers =    ['o',                's',        '^',      'D']

fig, ax = plt.subplots(figsize=(3.5, 3.0))
fig.subplots_adjust(left=0.15, right=0.85, top=0.88, bottom=0.12)

x = [0, 1]

for i, (name, s, io, c, m) in enumerate(zip(evaluators, setup, imageonly, colors, markers)):
    ax.plot(x, [s, io], '-', color=c, linewidth=1.5, zorder=3)
    ax.plot(x, [s, io], m, color=c, markersize=5, zorder=4)
    # Only percentages on sides — names go in legend
    ax.text(-0.05, s, f'{s:.0f}', ha='right', va='center', fontsize=6, color=c)
    ax.text(1.05, io, f'{io:.0f}', ha='left', va='center', fontsize=6, color=c)

ax.set_xlim(-0.25, 1.25)
ax.set_ylim(25, 95)
ax.set_xticks(x)
ax.set_xticklabels(['Setup-\nconditioned', 'Image-\nonly'], fontsize=7, ha='center')
ax.set_ylabel('Accuracy (%)', fontsize=7)
ax.tick_params(axis='y', labelsize=6)
ax.yaxis.grid(True, alpha=0.2, linewidth=0.5)
ax.set_axisbelow(True)

for spine in ['top', 'right', 'bottom']:
    ax.spines[spine].set_visible(False)
ax.spines['left'].set_alpha(0.3)
ax.tick_params(bottom=False)

# Legend with markers
handles = [mlines.Line2D([], [], color=c, marker=m, markersize=4, linewidth=1.2, label=n)
           for n, c, m in zip(evaluators, colors, markers)]
ax.legend(handles=handles, fontsize=5.5, loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=2, frameon=False, columnspacing=1.0, handletextpad=0.4)

fig.savefig(OUT / 'fig_result_rank_reversal.pdf', format='pdf')
fig.savefig(OUT / 'fig_result_rank_reversal.png', format='png', dpi=600)
plt.close()
print("Saved: fig_result_rank_reversal.pdf/.png")
