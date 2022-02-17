import figures
import matplotlib.pyplot as plt
from matplotlib import cm

# Circos plot
results = figures.load_results('significant')
G = figures.results_to_networkx(results)
fig = plt.figure(figsize=(10, 12))
ax1 = fig.add_subplot(121)
figures.plot_circos(G, ax=ax1)
ax2 = fig.add_subplot(122)
figures.plot_circos(G, ax=ax2, sex='male', title='Males')
fig.tight_layout()
fig.savefig('../Results/Plots/Figure2.pdf')

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
