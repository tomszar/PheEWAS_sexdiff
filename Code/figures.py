import networkx as nx
import pandas as pd
import nxviz as nv
import numpy as np
import matplotlib.pyplot as plt
from nxviz import annotate


def load_results(type: str = 'all'):
    '''
    Load results table and return only significant associations

    Parameters
    ----------
    type: str
        the type of results requested, it can be 'all', 'significant',
        'male', or 'female'

    Returns
    ----------
    results: pd.DataFrame
        results table
    '''
    results   = pd.read_csv('../Results/ResultsTable.csv')
    converged_columns = ['Converged_df',
                         'Converged_rf',
                         'Converged_dm',
                         'Converged_rm']
    converged = results.loc[:,converged_columns].all(axis=1)
    results = results.loc[converged]
    if type == 'significant':
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
            G.nodes[n]['name'] = (results[rows]['Variable_Category']).\
                                                              unique()[0]
            G.nodes[n]['order'] = order_by[G.nodes[n]['category']]
        else:
            G.nodes[n]['group'] = 'phenotype'
            ind = ph_descript[ph_descript[0] == n].index.to_list()[0]
            G.nodes[n]['category'] = ph_descript.iloc[ind,2]
            G.nodes[n]['name'] = ph_descript.iloc[ind,1]
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
                    pos[node] = np.array([x, y])
                    xd, yd = nv.polcart.to_cartesian(r=radius * 1.1,
                                    theta=nv.geometry.item_theta(nodes, node))
                    deg = nv.polcart.to_degrees(\
                                     nv.geometry.item_theta(nodes,
                                                            node))
                    degrees[node] = np.array([xd, yd, deg])

    else:
        for node, data in nt.iterrows():
            x, y = nv.polcart.to_cartesian(r=radius,
                                theta=nv.geometry.item_theta(nodes, node))
            pos[node] = np.array([x, y])
            xd, yd = nv.polcart.to_cartesian(r=radius * 1.1,
                                theta=nv.geometry.item_theta(nodes, node))
            deg = nv.polcart.to_degrees(\
                             nv.geometry.item_theta(nodes,
                                                    node))
            degrees[node] = np.array([xd, yd, deg])
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
            x = degrees[r[0]][0]
            y = degrees[r[0]][1]
            deg = degrees[r[0]][2]
            if -90 <= deg <= 90:
                rot = deg
            else:
                rot = deg - 180
            ha, va = annotate.text_alignment(x, y)
            ax.annotate(xy = (x, y),
                        ha = ha,
                        va = 'center',
                        text = r[1]['name'],
                        rotation=rot,
                        rotation_mode='anchor',
                        size=6)

    fontchanges = {'fontsize': 20}
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

def plot_miami(results,
               to_annotate: dict = None):
    '''
    Miami plot with female and male results.

    Parameters
    ----------
    results: pd.DataFrame
        results dataframe with male and female results
    to_annotate: dict
        dictionary with variable and outcome to annotate
    '''
    offset = 5  # Spacing between categories
    df = results.copy()
    
    # Transform pvalues
    df['log_pval_female'] = -np.log10(df['pvalue_female'])
    df['log_pval_male'] = np.log10(df['pvalue_male'])

    # Generating x position
    df['category_x_offset'] = df.groupby('Variable_Category').ngroup() * offset
    df['xpos'] = df.groupby(['Variable_Category', 'Variable_Name']).ngroup() + \
                 df['category_x_offset']
    
    # X labels
    x_labels = []
    x_labels_pos = []

    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    pval_columns = ['log_pval_female',
                    'log_pval_male']
    for category_num, (category_name, category_data) in enumerate(
            df.groupby('Variable_Category') ):
        if category_num % 2 == 0:
            base_color = 'lightgray'
        else:
            base_color = 'darkgray'
        left = category_data['xpos'].min() - (offset / 2) - 0.5
        right = category_data['xpos'].max() + (offset / 2) + 0.5
        width = right - left
        x_labels.append(category_name)
        x_labels_pos.append(left + width / 2)

        for pval_col in pval_columns:
            if pval_col == 'log_pval_female':
                ax = ax1
            else:
                ax = ax2
            ax.scatter(x='xpos',
                       y=pval_col,
                       label=None,
                       color=base_color,
                       s=30000 / len(df),
                       data=category_data,
                       alpha=0.6)
            # Annotate
            if to_annotate is not None:
                for var in to_annotate:
                    d = df.query('(Variable == @var) \
                                  and Outcome == @to_annotate[@var]')
                    x = d['xpos'].iloc[0]
                    y = d[pval_col].max()
                    name = d['Variable_Name'].iloc[0] + \
                           '\n' + \
                           d['Outcome_Name'].iloc[0]

                    if var == 'LBXCOT':
                        yoffset = 100
                    elif var == 'SMD070':
                        yoffset = 200
                    elif var == 'SMD090':
                        yoffset = 100
                    else:
                        yoffset = 150
                    
                    if pval_col == 'log_pval_female':
                        xytext = (x, y + yoffset)
                    else:
                        xytext = (x, y - yoffset)
                    ax.annotate(text=name,
                                xy=(x,y),
                                xytext=xytext,
                                ha='left',
                                arrowprops=dict(arrowstyle='->'))

            # Highlight sex differences
            for diff_num, (diff_name, diff_data) in enumerate(
                category_data.groupby('difference_type')):
                if diff_name == 'Pure':
                    high_color = 'orange'
                elif diff_name == 'Qualitative':
                    high_color = 'green'
                elif diff_name == 'Quantitative':
                    high_color = 'blue'
                
                if diff_name != 'None':
                    ax.scatter(x='xpos',
                               y=pval_col,
                               label=None,
                               color=high_color,
                               s=300000 / len(df),
                               data=diff_data,
                               alpha=0.7)

    labels = df['Variable_Category'].unique()
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        if ax == ax1:
            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels, 
                               rotation=45,
                               ha='right',
                               fontdict={'fontsize':8})
        else:
            ax.get_xaxis().set_visible(False)
        #ax.spines['left'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('../Results/Plots/Figure1.pdf')
    plt.show()

def plot_forest(results,
                colors):
    '''
    Create a forest (dot) plot showing the effect sizes
    between sexes

    Parameters
    ----------
    results: pd.DataFrame
        significant results to show
    colors: matplotlib.colors.ListedColormap
        colormap to use. Only takes the first two colors for female and male
        in that order
    '''
    fig = plt.figure(figsize=(8,12))
    ax = fig.add_subplot(111)
    sexes = ['female',
             'male']
    differences = ['Pure',
                   'Quantitative',
                   'Qualitative']
    df = results.sort_values(by=['difference_type',
                                 'Variable',
                                 'Outcome'])
    y = list(range(0, len(df)))
    beta_female = df.loc[:,'Beta_female']
    beta_male = df.loc[:,'Beta_male']
    ax.vlines([0], -1, 
              len(df) + 1,
              colors='k')
    ax.plot(beta_female,
            y,
            'bo')
    ax.plot(beta_male,
            y,
            'mo')
    
    plt.show()
    
