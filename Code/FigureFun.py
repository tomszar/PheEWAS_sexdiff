# Function for figures

def estimate_corr(dat):
    """
    Generate correlation and outputs a numpy array with zeros in the diagonal
    """
    import numpy as np
    import pandas as pd
    
    corr = dat.corr(method='spearman')
    corr = np.array(corr)
    np.fill_diagonal(corr, 0)
    return(corr)

def cluster_corr(corr_array, inplace=False, returnindex=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to each other 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    import scipy.cluster.hierarchy as sch
    import numpy as np
    import pandas as pd
    
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if returnindex:
        return idx
    else:
        if not inplace:
            corr_array = corr_array.copy()
        
        if isinstance(corr_array, pd.DataFrame):
            return corr_array.iloc[idx, :].T.iloc[idx, :]
        return corr_array[idx, :][:, idx]


def correl_plot(correlations, figs=(8,8), title=None, ticknames=None, filename=None, cm='coolwarm', vs=None):
    #Plot correlation overall
    import matplotlib.pyplot as plt
    import pandas as pd
    
    fig = plt.figure(figsize=figs)
    ax  = fig.add_subplot(111)
    im  = ax.imshow(pd.DataFrame(correlations), cmap=cm,)
        
    if ticknames is not None:
        plt.xticks(range(len(ticknames)), ticknames, rotation=90);
        plt.yticks(range(len(ticknames)), ticknames);
        
    fig.colorbar(im)
    ax.set_title(title)
    
    #Save figure
    if filename is not None:
        plt.savefig(filename)
    
    #Show
    if filename is None:
        plt.show()

def plot_betas(dat, exposure_name, beta_names=['Beta_f', 'Beta_m'], ci_names=['ci_f', 'ci_m'], vline_limits=[0, 10], title='', figs=(8,8)):
    '''
    Plot beta coefficients between two groups (e.g. males and females), defined by a particular exposure
    
    Parameters
    ----------
    - dat: DataFrame
        database where the beta values, and CI will be obtained
    - beta_names: list[str, str] (default ['Beta_f', 'Beta_m'])
        the column names of the beta coefficients to compare
    - ci_names: list of str (default ['ci_f', 'ci_m'])
        the column name of the intervals to use. Preferentially use confidence intervals
    - exposure_name: str
        the name of the exposure to use
    - vline_limits: list[int, int] (default [0, 10])
        the limits of the vertical line that passes through zero
    - title: str (default '')
        title name
    - figs: tuple(int, int) (default (8,8))
        figure size
        
    Returns
    -------
    None
    
    '''
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=figs, dpi=300)
    plt.errorbar(x = dat.loc[exposure_name, beta_names[0]], 
                 y = dat.loc[exposure_name].index, 
                 xerr = dat.loc[exposure_name, ci_names[0]],
                 fmt='.', label='Females')
    
    plt.errorbar(x = dat.loc[exposure_name, beta_names[1]], 
                 y = dat.loc[exposure_name].index, 
                 xerr = dat.loc[exposure_name, ci_names[1]], 
                 fmt='.', label='Males')
    
    plt.vlines(0, vline_limits[0], vline_limits[1], linestyle='dashed', linewidth=1, colors='black')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=-0.5)
    
    plt.title(title)
    plt.show()
    
def plot_effect_sizes(dat, pval_names=['pvalue_f', 'pvalue_m'], beta_names=['Beta_f', 'Beta_m'], ci_names=['ci_f', 'ci_m'], title='', figs=(8,8), pt=0.05,
                      x_limits=None):
    '''
    Plot effect sizes and p-values across a range of phenotype-exposures
    Parameters
    ----------
    - dat: DataFrame
        database where the beta values, and CI will be obtained
    - pval_names: list[str, str] (default ['pvalue_f', 'pvalue_m'])
        the column names of the pvalues to compare
    - beta_names: list[str, str] (default ['Beta_f', 'Beta_m'])
        the column names of the beta coefficients to compare
    - ci_names: list of str (default ['ci_f', 'ci_m'])
        the column name of the intervals to use. Preferentially use confidence intervals
    - title: str (default '')
        title name
    - figs: tuple(int, int) (default (8,8))
        figure size
    - pt: int (default 0.05)
        pvalue threshold to show
    - x_limits: tuple(int, int) (default None)
        set limits to x axis in second plot
        
    Returns
    -------
    None
    
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    
    fig = plt.figure(figsize=figs)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    old_e = []
    pos   = 0
    step  = 1
    lowylim = -1
    ax1.set_ylim(lowylim, len(dat))
    ax2.set_ylim(lowylim, len(dat))
    ax3.set_ylim(lowylim, len(dat))
    if x_limits is not None:
        ax2.set_xlim(x_limits)
    
    for i in range(0, len(dat)):
        new_e = dat.sort_index().index[i][0]
        new_p = dat.sort_index().index[i][1]
        if new_e != old_e:
            ax1.text(0, pos, new_e)
        ax1.text(0.5, pos, new_p)
        old_e = new_e
        pos   = pos + step
        
    ax2.errorbar(x = dat.sort_index().loc[:, beta_names[0]],
                 y = list(range(0,len(dat))), 
                 xerr = dat.sort_index().loc[:, ci_names[0]],
                 fmt='o', label='Females')
    
    ax2.errorbar(x = dat.sort_index().loc[:, beta_names[1]], 
                 y = list(range(0,len(dat))), 
                 xerr = dat.sort_index().loc[:, ci_names[1]],
                 fmt='o', label='Males')

    ax2.vlines(0, 0, len(dat), linestyle='dashed', linewidth=1, colors='black')
    
    ax3.plot(-np.log10(np.array(dat.sort_index().loc[:, pval_names[0]])),
             list(range(0,len(dat))), 'o')
    
    ax3.plot(-np.log10(np.array(dat.sort_index().loc[:, pval_names[1]])), 
             list(range(0,len(dat))), 'o')
    
    ax3.vlines(0, 0, len(dat), linestyle='dashed', linewidth=1, colors='black')
    ax3.vlines(0, -np.log10(pt), len(dat), linestyle='dashed', linewidth=1, colors='red')
    
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')
    
    ax2.set_title('Beta coefficients')
    ax3.set_title('P-values')
    
    ax2.legend(bbox_to_anchor=(1.8, 0.95), loc='upper left', borderaxespad=-0.5)
    fig.suptitle(title, y=0.92)
    plt.show()
    
#def correl_two_plots(dfs, figs=(8,8), titles=None, ticknames=None, filenames=None):
    