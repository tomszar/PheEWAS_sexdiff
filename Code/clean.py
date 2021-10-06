import os
import clarite
import pandas as pd

class NhanesRaw:
    '''
    NHANES raw dataset class
    '''
    def __init__(self):
        '''
        Initiate the class

        Attributes
        ----------
        data: pd.DataFrame
            NHANES raw data with only adults (>18yo)
        var_description: dict
            Dictionary with the variables description
        var_category: dict
            Dictionary with the variables category
        phenotypes: list of str
            List of phenotype variable names
        exposures: list of str
            List of exposure variables names
        covariates: list of str
            List of covariate variables names
        cohorts: list of str
            List of cohort names
        cohorts_bool: list of bool
            List of boolean with the cohort that each participant belongs
        '''
        print('-----Loading NHANES raw data-----')

        # Setting data path
        data_path = '../Data/'

        # Loading data
        var_info  = _load_var_information(data_path)
        self.var_description = var_info[0]
        self.var_category    = var_info[1]

        filename = os.path.join(data_path,
                                'nh_99-06',
                                'MainTable.csv')
        nhanes   = pd.read_csv(filename).\
                      rename(columns={'SEQN':'ID'}).\
                      set_index('ID')
        self.data = nhanes

        # Dividing data
        self.covariates = ['black',
                           'mexican',
                           'other_hispanic',
                           'other_eth',
                           'SES_LEVEL',
                           'RIDAGEYR',
                           'SDDSRVYR',
                           'BMXBMI',
                           'female',
                           'male']
        
        exposure_categories = ['alcohol use',
                               'bacterial infection',
                               'cotinine',
                               'diakyl',
                               'dioxins',
                               'food component recall',
                               'furans',
                               'heavy metals',
                               'housing',
                               'hydrocarbons',
                               'nutrients',
                               'occupation',
                               'pcbs',
                               'perchlorate',
                               'pesticides',
                               'phenols',
                               'phthalates',
                               'phytoestrogens',
                               'polybrominated ethers',
                               'polyflourochemicals',
                               'sexual behavior',
                               'smoking behavior',
                               'smoking family',
                               'social support',
                               'street drug',
                               'sun exposure',
                               'supplement use',
                               'viral infection',
                               'volatile compounds']
        
        phenotype_categories = ['biochemistry',
                                'blood',
                                'hormone']
        
        self.phenotypes = self._get_var_from_category(phenotype_categories)
        self.exposures  = self._get_var_from_category(exposure_categories)

        # Selecting participants
        # Only adults
        keep_adults = self.data['RIDAGEYR'] >= 18
        self.data   = self.data.loc[keep_adults]
        # Without missing covariates
        self.data = clarite.modify.\
                            rowfilter_incomplete_obs(self.data,
                                                     only=self.covariates)

        # Creating cohorts
        discovery   = (self.data['SDDSRVYR']==1) |\
                      (self.data['SDDSRVYR']==2) 
        replication = (self.data['SDDSRVYR']==3) |\
                      (self.data['SDDSRVYR']==4)
        
        females = self.data['female'] == 1
        males   = self.data['male'] == 1

        self.cohorts = ['discovery females',
                        'discovery males',
                        'replication females',
                        'replication males']
        
        self.cohorts_bool = [discovery & females,
                             discovery & males,
                             replication & females,
                             replication &  males]

        print('Loaded NHANES data')
        self._print_summary()
        print('')

    def _print_summary(self):
        '''
        Print a summary information with number of participants in each cohort,
        phenotypes, exposures, and covariates
        '''
        n_cov   = len(self.covariates)
        n_pheno = len(self.phenotypes)
        n_expo  = len(self.exposures)
        n_part  = len(self.data)

        print('There are ' +
              str(n_part) + 
              ' participants, ' +
              str(n_cov) + 
              ' covariates, ' +
              str(n_pheno) + 
              ' phenotypes, and ' +
              str(n_expo) + 
              ' exposures')

        print('For each cohort, there are the following sample sizes:')
        for i,l in enumerate(self.cohorts):
            n = len(self.data[self.cohorts_bool[i]])
            print(l + 
                  ': ' + 
                  str(n))

    def _get_var_from_category(self,
                               var_categories:list):
        '''
        From a list of variable categories, retrieve a list
        of variable names
    
        Parameters
        ----------
        var_categories: list of str
            category names
    
        Returns
        ----------
        var_names: list of str
            variable names
        '''
        var_names = []
        for i in range(len(var_categories)):
            keys_list = [k for k, 
                         v in self.var_category.items()\
                         if v == var_categories[i]]
            for l in keys_list:
                var_names.append(l)
        return(var_names)

def _load_var_information(datapath:str):
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

