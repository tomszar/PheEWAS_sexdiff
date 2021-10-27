import os
import pandas as pd

def load_var_information(datapath:str='../Data/'):
    '''
    Load variable description file in a suitable dictionary,
    with either the variable description or the variable category

    Parameters
    ----------
    datapath: str
        folder path for the data location

    Returns
    ----------
    var_info: list of dict
        variable information with description and category
    '''
    filename = os.path.join(datapath, 'nh_99-06', 'VarDescription.csv')
    var_description = pd.read_csv(filename)\
                        .drop_duplicates()\
                        .set_index('var')

    # Convert variable descriptions to a dictionary for convenience
    var_description_dict = var_description['var_desc'].to_dict()
    var_category_dict    = var_description['category'].to_dict()

    # Output list of dict
    var_info = [var_description_dict,
                var_category_dict]
    return(var_info)

def get_covariate_list():
    '''
    Get the list of covariates

    Returns
    ----------
    covariates: list of str
        list of covariate names
    '''
    covariates = ['black',
                  'mexican',
                  'other_hispanic',
                  'other_eth',
                  'SES_LEVEL',
                  'RIDAGEYR',
                  'SDDSRVYR',
                  'BMXBMI',
                  'female',
                  'male']
    return(covariates)

def get_cohorts(nhanes_data):
        '''
        Get cohorts based on survey year and sex from nhanes dataset

        Parameters
        ----------
        nhanes_data: pd.DataFrame
            NHANES dataframe with survey year and sex information

        Returns
        ----------
        cohorts: list of str
            List of cohort names
        cohorts_bool: list of bool
            List of boolean with the cohort that each participant belongs
        '''
        discovery   = (nhanes_data['SDDSRVYR']==1) |\
                      (nhanes_data['SDDSRVYR']==2) 
        replication = (nhanes_data['SDDSRVYR']==3) |\
                      (nhanes_data['SDDSRVYR']==4)
        
        females = nhanes_data['female'] == 1
        males   = nhanes_data['male'] == 1

        cohorts = ['discovery females',
                   'discovery males',
                   'replication females',
                   'replication males']
        
        cohorts_bool = [discovery & females,
                        discovery & males,
                        replication & females,
                        replication &  males]

        return([cohorts,
                cohorts_bool])

def add_var_from_dict(index_var, 
                      var_dict):
    '''
    Add new variable from dict

    Parameters
    ----------
    index_var
    '''
    res = []
    for i in range(len(index_var)):
        if index_var[i] in var_dict.keys():
            res.append(var_dict[index_var[i]])
        else:
            res.append('')
    return(res)
