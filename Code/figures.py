import networkx as nx
import pandas as pd
import nxviz as nv
import numpy as np
import matplotlib.pyplot as plt
from nxviz import annotate


def load_results():
    '''
    Load results table and return only significant associations

    Returns
    ----------
    results: pd.DataFrame
        results table
    '''
    results   = pd.read_csv('../Results/ResultsTable.csv')
    keep_rows = ~(results['difference_type'] == 'None')
    results   = results[keep_rows]
    return results

def results_to_networkx(results):
    '''
    Transforms the results table into a network graph with
    edge_values and group

    Parameters
    ----------
    results: pd.DataFrame
        results table

    Returns
    ----------
    G: nx.Graph
        network graph 
    '''
    order_by = {'cotinine': 0,
                'smoking behavior': 1,
                'smoking family': 2,
                'heavy metals': 3,
                'pcbs': 4, 
                'volatile compounds': 5,
                'phthalates': 6,
                'dioxins': 7,
                'food component recall': 8,
                'supplement use': 9,
                'nutrients': 10,
                'albumin': 11,
                'AP': 12,
                'GGT': 13,
                'electrolyte': 14,
                'glucose': 15,
                'homocysteine': 16,
                'inflammation': 17,
                'iron': 18,
                'kidney function': 19,
                'lipids': 20,
                'osmolality': 21,
                'red blood cells': 22,
                'total protein': 23,
                'MMA': 24,
                'white blood cells': 25}

    ph_descript = pd.read_csv('../Data/phenotype_descriptions.csv',
                              header=None)

    G = nx.from_pandas_edgelist(results, 
                                source = 'Variable', 
                                target = 'Outcome')
    # Whether the association is significant
    threshold = 0.05
    for v, o in G.edges:
        out = results['Outcome']  == o# o and v are outcome and variables
        var = results['Variable'] == v
        row = results[var & out]
        if len(row) == 0:
            out = results['Outcome']  == v# o and v are outcome and variables
            var = results['Variable'] == o
            row = results[var & out]

        if float(row['pvalue_female']) < threshold:
            G.edges[v, o]['significant_female'] = 1
        else:
            G.edges[v, o]['significant_female'] = 0
        
        if float(row['pvalue_male']) < threshold:
            G.edges[v, o]['significant_male'] = 1
        else:
            G.edges[v, o]['significant_male'] = 0
        
        G.edges[v, o]['beta_female'] = float(row['Beta_female'])
        G.edges[v, o]['beta_male'] = float(row['Beta_male'])
    
    for n in G.nodes:
        if n in list(results['Variable']):
            rows = results['Variable'] == n
            G.nodes[n]['group'] = 'exposure'
            G.nodes[n]['category'] = (results[rows]['Variable_Category']).\
                                                              unique()[0]
            G.nodes[n]['order'] = order_by[G.nodes[n]['category']]
        else:
            G.nodes[n]['group'] = 'phenotype'
            ind = ph_descript[ph_descript[0] == n].index.to_list()[0]
            G.nodes[n]['category'] = ph_descript.iloc[ind,2]
            G.nodes[n]['order'] = order_by[G.nodes[n]['category']]

    return G

def novel_circos(nt: pd.DataFrame,
                 group_by = None,
                 sort_by = None,
                 radius: float = None):
    '''
    Circos layout with separations between groups defined in group_by
    '''
    pos = dict()
    degrees = dict()
    nt = nv.utils.group_and_sort(nt,
                                 group_by,
                                 sort_by)
    nodes = list(nt.index)

    if radius is None:
        radius = nv.geometry.circos_radius(len(nodes))
    
    if group_by:
        # Add empty elements at the end of each group to separate them
        df2 = pd.DataFrame.from_dict({'category': ['empty'],
                                      'group': ['empty'],
                                      'size': [0]})
        df2 = pd.concat([df2]*3,
                        ignore_index=True)
        xt = pd.DataFrame(columns = nt.columns)    
        for grp, df in nt.groupby(group_by):
            df = df.append(df2)
            xt = xt.append(df)

        nodes = list(xt.index)
        radius = nv.geometry.circos_radius(len(nodes))

        for grp, df in xt.groupby(group_by):
            if grp != 'empty':
                for node, data in df.iterrows():
                    x, y = nv.polcart.to_cartesian(r=radius,
                                    theta=nv.geometry.item_theta(nodes, node))
                    degrees[node] = nv.polcart.to_degrees(\
                                       nv.geometry.item_theta(nodes,
                                                              node))
                    pos[node] = np.array([x, y])

    else:
        for node, data in nt.iterrows():
            x, y = nv.polcart.to_cartesian(r=radius,
                                theta=nv.geometry.item_theta(nodes, node))
            degrees[node] = nv.polcart.to_degrees(nv.geometry.item_theta(nodes,
                                                                         node))
            pos[node] = np.array([x, y])
    return pos, degrees

def plot_circos(G,
                ax,
                title='Females',
                sex='female'):
    '''
    Generate circos plot
    '''
    nt = nv.utils.node_table(G)
    nt['size'] = 1
    et = nv.utils.edge_table(G)
    
    ax = plt.gca()
    pos, degrees = novel_circos(nt,
                                group_by='group',
                                sort_by='order')
    
    node_color = group_colormap(nt['category'])
    alpha = nv.nodes.transparency(nt, alpha_by=None)
    size = nv.nodes.node_size(nt, 'size')
    patches = nv.nodes.node_glyphs(nt,
                                   pos,
                                   node_color=node_color,
                                   alpha=alpha,
                                   size=size)
    for patch in patches:
        ax.add_patch(patch)
    
    # Customize edge styling
    edge_color = nv.edges.edge_colors(et,
                                      nt=None,
                                      color_by=None,
                                      node_color_by=None)
    beta = 'beta_' + sex
    significant = 'significant_' + sex
    edge_color[et[beta] < 0] = 'red'
    edge_color[et[beta] > 0] = 'blue'
    lw = np.sqrt(np.abs(et[beta])) * 8
    lw[et[significant] == 0] = 0
    alpha = nv.edges.transparency(et, alpha_by=None)
    patches = nv.lines.circos(et,
                              pos,
                              edge_color=edge_color,
                              alpha=alpha,
                              lw=lw,
                              aes_kw={'fc': 'none'})
    for patch in patches:
        ax.add_patch(patch)
    
    # Annotation for categories in exposures,
    # and names for phenotypes
    for grp, df in nt.groupby('category'):
        if df['group'][0] != 'phenotype':
            ind = list(df.index)
            xys = []
            for i in ind:
                val = tuple(pos.get(i))
                xys.append(val)
            
            ind = round(len(xys) / 2)
            if ind == 1:
                ind = 0
        
            if grp == 'phthalates':
                offsety = 1.4
            else:
                offsety = 0
            x = ( (xys[ind][0]) * 1.05 )
            y = ( (xys[ind][1]) * 1.05 ) + offsety
            ha, va = annotate.text_alignment(x, y)
            ax.annotate(grp, xy=(x, y), ha=ha, va=va)
    
    for r in nt.iterrows():
        if r[1]['group'] == 'phenotype':
            ax.text(x = pos[r[0]][0],
                    y = pos[r[0]][1],
                    s = r[0],
                    rotation = degrees[r[0]])

    fontchanges = {'fontsize': 20}
                   #'fontweight': 50}

    ax.set_title(title,
                 fontdict=fontchanges,
                 pad=25)
    
    nv.plots.rescale(G)
    nv.plots.aspect_equal()
    nv.plots.despine()

def group_colormap(data: pd.Series):
    a = plt.get_cmap('tab20')
    b = plt.get_cmap('tab20b')
    c = plt.get_cmap('tab20c')
    cmap = {'cotinine': a(0),
            'smoking behavior': a(1),
            'smoking family': a(2),
            'heavy metals': a(3),
            'pcbs': a(4), 
            'volatile compounds': a(5),
            'phthalates': a(6),
            'dioxins': a(7),
            'food component recall': a(8),
            'supplement use': a(9),
            'nutrients': a(10),
            'albumin': b(0),
            'AP': b(1),
            'GGT': b(2),
            'electrolyte': b(3),
            'glucose': c(0),
            'homocysteine': c(1),
            'inflammation': c(2),
            'iron': c(3),
            'kidney function': c(12),
            'lipids': c(13),
            'osmolality': c(14),
            'red blood cells': c(15),
            'total protein': b(16),
            'MMA': b(17),
            'white blood cells': b(18)}
    return data.apply(lambda x: cmap.get(x))
