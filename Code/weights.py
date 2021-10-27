import os
import pandas as pd

def load_weights_raw(datapath:str='../Data'):
    '''
    Load survey weight information,
    from the discovery and replication dataset.

    Parameters
    ----------
    datapath: str
        path to the weight data

    Returns
    ----------
    weights: list
        list of weight data
    '''
    files = ['weights_discovery.txt',
             'weights_replication.txt',
             'VarWeights.csv']
    weights = []
    for f in files:
        filename = os.path.join(datapath,
                                f)
        if f == 'VarWeights.csv':
            dat = pd.read_csv(filename).\
                     set_index('variable_name')
        else:
            dat = pd.read_csv(filename,
                              sep='\t').\
                                rename(columns={'SEQN':'ID'}).\
                                set_index('ID')
        weights.append(dat)
 
    # Divide the survey weights replication by 2 to get 4 year weights
    weights[1].iloc[:,3:] = \
        weights[1].iloc[:,3:] / 2
         
    weights[2] = _harmonize_weights_to_map(weights)
    return([weights[0],
            weights[1],
            weights[2]])

def load_weights_cleaned(datapath:str='../Results/'):
    '''
    Load cleaned weight information
    
    Parameters
    ----------
    datapath: str
        path to folder where the weight data is stored
    
    Returns
    ----------
    weights: list
        list of dataframes of discovery, replication, and map data
    '''
    files = ['weights_discovery.csv',
             'weights_replication.csv',
             'weights_maps.csv']

    weights = []

    for f in files:
        filename = os.path.join(datapath,
                                f)
        if f == 'weights_maps.csv':
            idx_name = 'variable_name'
        else:
            idx_name = 'ID'
        dat = pd.read_csv(filename).\
                 set_index(idx_name)
        weights.append(dat)
    
    return(weights)

def _harmonize_weights_to_map(weights:list):
        '''
        Remove the weights variables that don't coincide between map and data

        Parameters
        ----------
        mapped: pd.DataFrame
            Updated map data
        '''
        print('-----Internal harmonizing of survey weights-----')
        cohorts = ['discovery',
                   'replication']
        remove_vars = []
        mapped = weights[2]
        for i,c in enumerate(cohorts):
            cols = weights[i]
            vals = set(mapped[cohorts[i]])

            # Which survey weights are not present
            weights_absent = []
            for l in vals:
                if l not in cols:
                    weights_absent.append(l)

            # Which variables need those survey weights
            for key, value in enumerate(mapped[cohorts[i]]):
                for m in weights_absent:
                    if m == value:
                        remove_vars.append(mapped[cohorts[i]].index[key])
           
        mapped = mapped.drop(remove_vars)
        return(mapped)
