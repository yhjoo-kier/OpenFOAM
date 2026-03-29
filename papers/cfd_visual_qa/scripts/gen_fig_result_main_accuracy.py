"""Fig 4: (a) Accuracy bars 4 models x 2 conditions + (b) setup dependency."""
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

# 3-trial data
claude_s = [90.0, 90.0, 86.7]; claude_io = [33.3, 33.3, 33.3]
gpt_s = [56.7, 70.0, 76.7]; gpt_io = [43.3, 56.7, 46.7]
gemini_s = [56.7, 53.3, 63.3]; gemini_io = [63.3, 56.7, 40.0]
expert_s = [73.8]; expert_io = [66.7]

names = ['Claude', 'GPT-5.4', 'Expert', 'Gemini']
setup_mean = [np.mean(claude_s), np.mean(gpt_s), np.mean(expert_s), np.mean(gemini_s)]
io_mean = [np.mean(claude_io), np.mean(gpt_io), np.mean(expert_io), np.mean(gemini_io)]
setup_std = [np.std(claude_s), np.std(gpt_s), 0, np.std(gemini_s)]
io_std = [np.std(claude_io), np.std(gpt_io), 0, np.std(gemini_io)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.6),
    gridspec_kw={'width_ratios': [1.2, 0.8], 'wspace': 0.40})

# ── (a) Grouped bar chart — uniform colors for conditions ──
x = np.arange(len(names))
w = 0.30
c_setup = '#2B5EA7'
c_io = '#B0C8E8'

bars1 = ax1.bar(x - w/2, setup_mean, w, yerr=setup_std, capsize=2,
                color=c_setup, edgecolor='white', linewidth=0.5, label='Setup-conditioned')
bars2 = ax1.bar(x + w/2, io_mean, w, yerr=io_std, capsize=2,
                color=c_io, edgecolor='white', linewidth=0.5, label='Image-only')

for bar, std in zip(bars1, setup_std):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 1.5,
             f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=5)
for bar, std in zip(bars2, io_std):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 1.5,
             f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=5)

ax1.set_ylabel('Accuracy (%)', fontsize=7)
ax1.set_ylim(0, 102)
ax1.set_xticks(x)
ax1.set_xticklabels(names, fontsize=6)
ax1.tick_params(labelsize=6)
ax1.legend(fontsize=5.5, loc='upper right', frameon=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.text(-0.14, 1.03, '(a)', transform=ax1.transAxes, fontsize=7, fontweight='semibold')

# ── (b) Setup dependency horizontal bars ──
delta = [s - io for s, io in zip(setup_mean, io_mean)]
sorted_idx = np.argsort(delta)[::-1]

y_pos = np.arange(len(names))
ax2.barh(y_pos, [delta[i] for i in sorted_idx],
         color='#888888', edgecolor='white', height=0.55)

for j, i in enumerate(sorted_idx):
    val = delta[i]
    sign = '+' if val > 0 else ''
    ax2.text(val + 1.5, j, f'{sign}{val:.1f}', ha='left', va='center', fontsize=5.5)

ax2.set_yticks(y_pos)
ax2.set_yticklabels([names[i] for i in sorted_idx], fontsize=6)
ax2.set_xlabel('\u0394 accuracy (pp)', fontsize=7)
ax2.axvline(0, color='black', linewidth=0.5, alpha=0.3)
ax2.tick_params(labelsize=6)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.text(-0.22, 1.03, '(b)', transform=ax2.transAxes, fontsize=7, fontweight='semibold')

fig.savefig(OUT / 'fig_result_main_accuracy.pdf', format='pdf', bbox_inches='tight')
fig.savefig(OUT / 'fig_result_main_accuracy.png', format='png', dpi=600, bbox_inches='tight')
plt.close()
print("Saved: fig_result_main_accuracy.pdf/.png")
