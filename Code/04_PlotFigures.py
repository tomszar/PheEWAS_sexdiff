import figures
import nxviz as nv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm


# Circos plot
types = ['Pure',
         'Quantitative',
         'Qualitative',
         'significant']
offsets = [None,
           {'food component recall': (0, 1.2),
            'heavy metals': (0, -1),
            'phthalates': (0, -0.5)},
           None,
           {'phthalates': (0, 1.2)}]

fig = plt.figure(figsize=(10, 16))
axs = gridspec.GridSpec(9, 5, figure=fig)
for i in range(4):
    results = figures.load_results(types[i])
    G = figures.results_to_networkx(results)
    first = (i*2) + 1
    last = ((i+1) * 2) + 1
    ax1 = fig.add_subplot(axs[first:last, 1:3])
    ax2 = fig.add_subplot(axs[first:last, 3:])
    ax_text = fig.add_subplot(axs[first:last, 0])
    figures.plot_circos(G, ax=ax1,
                        offsets=offsets[i])
    figures.plot_circos(G, ax=ax2,
                        sex='male',
                        offsets=offsets[i])
    if types[i] == 'significant':
        text = 'Total'
    else:
        text = types[i]
    stext = ax_text.text(1, 0.5, text,
                         size=18, rotation=90,
                         va='center')
    nv.plots.despine(ax=ax_text)

for sex in ['Females', 'Males']:
    if sex == 'Females':
        ax_sex = fig.add_subplot(axs[0, 1:3])
    elif sex == 'Males':
        ax_sex = fig.add_subplot(axs[0, 3:])
    ftext = ax_sex.text(0.5, -0.5, sex,
                        ha='center', size=18)
    nv.plots.despine(ax=ax_sex)

axs.tight_layout(fig, pad=0,
                 rect=(0.02, 0, 0.98, 1))
fig.savefig('../Results/Plots/Figure2.pdf',
            dpi=600)

# Dot plot
results = figures.load_results('significant')
color_pastel = cm.get_cmap('Set2')
figures.plot_forest(results,
                    color_pastel)

# Miami plot
results = figures.load_results()
figures.plot_miami(results,
                   to_annotate={'SMD070': 'LBXHCY',
                                'LBXCOT': 'LBXSAL',
                                'SMD090': 'LBXHCY',
                                'LBXRBF': 'LBXSAL'})
