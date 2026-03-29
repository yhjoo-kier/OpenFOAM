"""Fig 3: API-isolated evaluation protocol diagram."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
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

fig, ax = plt.subplots(figsize=(7.0, 2.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 2.2)
ax.set_axis_off()

# Stage box positions and labels
stages = [
    (0.3, 'CFD Cases\n60 simulations\n10 scenarios'),
    (2.3, 'Rendering\nstandardized\ncontour plots'),
    (4.3, 'Prompt\nsetup text\n+ image'),
    (6.3, 'API Call\nbase64 image\nno filesystem'),
    (8.3, 'Scoring\nOK / Anomaly\nvs ground truth'),
]

box_w, box_h = 1.6, 1.2
y_center = 1.3

for x, label in stages:
    box = FancyBboxPatch((x, y_center - box_h/2), box_w, box_h,
                          boxstyle='round,pad=0.08', facecolor='#EBF4FA',
                          edgecolor='#4A7FC1', linewidth=0.8)
    ax.add_patch(box)
    ax.text(x + box_w/2, y_center, label, ha='center', va='center',
            fontsize=6, color='#333333', linespacing=1.3)

# Arrows between stages
arrow_style = 'simple,head_width=4,head_length=3'
for i in range(len(stages) - 1):
    x1 = stages[i][0] + box_w
    x2 = stages[i+1][0]
    ax.annotate('', xy=(x2, y_center), xytext=(x1, y_center),
                arrowprops=dict(arrowstyle='->', color='#4A7FC1', lw=1.2))

# "x3 trials" annotation above API box
ax.text(6.3 + box_w/2, y_center + box_h/2 + 0.15, '\u00d73 independent trials',
        ha='center', va='bottom', fontsize=5.5, color='#2E7D32', style='italic')

# Bottom annotation: isolation barrier
ax.plot([5.9, 5.9], [y_center - box_h/2 - 0.3, y_center + box_h/2 + 0.05],
        color='#C62828', linewidth=1.5, linestyle='--')
ax.text(5.9, y_center - box_h/2 - 0.35, 'isolation\nbarrier', ha='center', va='top',
        fontsize=5, color='#C62828')

fig.savefig(OUT / 'fig_eval_protocol.pdf', format='pdf', bbox_inches='tight')
fig.savefig(OUT / 'fig_eval_protocol.png', format='png', dpi=600, bbox_inches='tight')
plt.close()
print("Saved: fig_eval_protocol.pdf/.png")
