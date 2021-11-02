import os
import utils
import weights
import clarite
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats

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
        weights_discovery: pd.DataFrame
            Data with weight information for the discovery
        weights_replication: pd.DataFrame
            Data with weight information for the replication
        weights_maps: pd.DataFrame
            Mapping data for variable - weight
        '''
        print('-----Loading NHANES raw data-----')

        # Setting data path
        data_path = '../Data/'

        # Loading data
        var_info  = utils.load_var_information(data_path)
        self.var_description = var_info[0]
        self.var_category    = var_info[1]

        filename = os.path.join(data_path,
                                'nh_99-06',
                                'MainTable.csv')
        nhanes   = pd.read_csv(filename).\
                      rename(columns={'SEQN':'ID'}).\
                      set_index('ID')
        self.data = nhanes

        # Adjusting lipids for statins
        self._adjust_lipids()

        # Dividing data
        self.covariates = utils.get_covariate_list()
        
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

        # Keeping only phenotypes, exposures, and covariates in data
        self._keep_data_from_lists()

        # Selecting participants
        # Only adults
        keep_adults = self.data['RIDAGEYR'] >= 18
        self.data   = self.data.loc[keep_adults]
        # Without missing covariates
        self.data = clarite.modify.\
                            rowfilter_incomplete_obs(self.data,
                                                     only=self.covariates)
        
        self._get_cohorts()

        # Weight data
        weight = weights.load_weights_raw('../Data')
        self.weights_discovery   = weight[0]
        self.weights_replication = weight[1]
        self.weights_maps        = weight[2]
        self._harmonize_weights_to_data()

        print('Loaded NHANES data')
        self._print_summary()
        print('')

    def categorize_variables(self):
        '''
        Categorize variables using clarite and remove constant variables,
        except covariates, done for each cohort independently.

        Returns
        ----------
        data: pd.DataFrame
            raw data with variables categorized and constant variables removed
        '''
        print('-----Categorizing data-----')
        
        # Make categorization in each and decide on how to achieve agreement
        to_remove = []
        var_types = []
        for i,b in enumerate(self.cohorts_bool):
            print('Categorizing in ' + 
                  self.cohorts[i] + 
                  ' cohort')
            temp_data = clarite.modify.categorize(self.data[b],
                                                  cont_min=10)
            var_types    = clarite.describe.get_types(temp_data)
            var_constant = var_types[var_types == 'constant'].index
            var_unknown  = var_types[var_types == 'unknown'].index
            to_remove.extend(var_constant)
            to_remove.extend(var_unknown)
            for c in temp_data.columns:
                if c not in self.data:
                    to_remove.append(c)
        
        to_remove = set(to_remove)
        print('Removing ' + 
              str(len(to_remove)) +
              ' variables (constant, only NAs ' + 
              'and continuous with less than 10 values)')
        self._remove_vars_in_lists(list(to_remove))
        
        self.data = clarite.modify.categorize(self.data,
                                              cont_min=10)
        self._update_lists_from_data()

        var_types    = clarite.describe.get_types(self.data)
        var_unknown  = var_types[var_types == 'unknown'].index 

        if len(var_unknown) > 0:
            print('Categorizing ' +
                  str(len(var_unknown)) +
                  ' unknown variables as continuous')
    
            self.data = clarite.modify.\
                                make_continuous(self.data,
                                                only=var_unknown)
        print('')        

    def clarite_filters(self):
        '''
        Run clarite filters in each cohort:
        - colfilter_min_n to remove variables with less than 200 non-Na values
        - colfilter_min_cat_n to remove variables with less than 200 values in
          a category
        - colfilter_percent_zero to remove variables with 90% of non-NA values
          equal to zero

        Returns
        ----------
        data: List of pd.Dataframe
            data with clarite filters applied
        '''
        print('-----Running CLARITE filters-----')
        self._get_cohorts()
        funcs = [clarite.modify.colfilter_min_n,
                 clarite.modify.colfilter_min_cat_n,
                 clarite.modify.colfilter_percent_zero]
        
        for fu in funcs:
            remove_vars = []
            for i in range(len(self.cohorts_bool)):
                print('Running in ' + 
                      self.cohorts[i] + 
                      ' cohort')
                temp_data = fu(self.data[self.cohorts_bool[i]],
                               skip=self.covariates)
                removed   = self._get_removed_vars(temp_data)
                remove_vars = remove_vars + list(removed)
            
            remove_vars = set(remove_vars)
            print('Removing ' +
                  str(len(remove_vars)) + 
                  ' variables across cohorts')
            self._remove_vars_in_lists(remove_vars)
            print('')
        self._print_summary()
        print('')

    def transform_continuous_log2(self):
        '''
        Log2 transform all continuous variables, and reflect those with
        negative skews in each cohort independently

        Returns
        ----------
        data: List of pd.Dataframe
            data with log2 transformed continuous variables
        '''
        self._get_cohorts()
        
        print('-----Log2 transformation-----')
        var_continuous  = self._get_data_type()
        print('Transforming ' + 
              str(len(var_continuous)) + 
              ' variables across cohorts')

        for i in range(len(self.cohorts_bool)):
            bools = self.cohorts_bool[i]
            
            for var in var_continuous:
                vmin = self.data.loc[bools, var].min()
                if vmin == 0: # add a small constant to zero values
                    dat_temp = np.array(self.data.loc[bools, var])
                    min_observed = np.nanmin(dat_temp[\
                                   np.nonzero(dat_temp)])
                    self.data.loc[bools, var] = self.data.loc[bools, var].\
                                                     replace(to_replace = 0,
                                                             value =
                                                             min_observed / 2)
                self.data.loc[bools, var] = \
                     np.log2(self.data.loc[bools, var])
        print('')

    def zscore_normalization(self):
        '''
        Apply z-score normalization to values to mean center and unit variance,
        by substracting whole column by mean, and dividing by the 
        standard deviation. 
        Done in each cohort independently

        Returns
        ----------
        data: List of pd.Dataframe
            data with scaled continuous variables
        '''
        self._get_cohorts
        
        print('-----Z-score normalization-----')
        var_continuous  = self._get_data_type()
        print('Scaling ' + 
              str(len(var_continuous)) + 
              ' variables across cohorts')

        for i in range(len(self.cohorts_bool)):
            bools = self.cohorts_bool[i]
            self.data.loc[bools, var_continuous] = \
                 self.data.loc[bools, var_continuous].\
                      apply(stats.zscore,
                            nan_policy='omit')

    def remove_continuous_outliers(self):
        '''
        Remove outliers for all continuous variables.
        If zscore normalization was done in each cohort 
        independently, removal of outliers will be as if was
        done in each cohort independently as well.

        Returns
        ----------
        data: list of pd.Dataframe
            data with continuous outliers removed
        '''
        print('-----Removing continuous outliers-----')
        var_continuous  = self._get_data_type()
        for c in self.covariates:
            if c in var_continuous:
                var_continuous.remove(c)

        self.data = clarite.modify.\
                           remove_outliers(self.data,
                                           only = var_continuous)

    def remove_perfect_correlations(self):
        '''
        Remove one random variable from a pair that have
        a perfect correlation
        '''
        print('-----Identifying identical variables-----')
        phenos_expos = self.phenotypes + \
                       self.exposures

        mat = np.array(self.data[phenos_expos].corr())

        #Replace diagonal
        np.fill_diagonal(mat, 0)
        
        #Identify variable
        bool_mat    = np.logical_or(mat == -1, mat == 1)
        where_index = np.where(bool_mat)

        #Remove vars
        remove_vars = []
        for i in range(0, len(where_index[0]), 2):
            var = phenos_expos[where_index[0][i]]
            remove_vars.append(var)

        print('Removing ' + 
              str(remove_vars))
        self._remove_vars_in_lists(remove_vars)
        print('')

    def plot_variables(self,
                       plotpath: str = '../Results/Plots/Inspect',
                       suffix:str = '',
                       only_continuous:bool = False):
        '''
        Plot variables

        Parameters
        ----------
        plotpath: str
            folder path to save the plots
        suffix: str
            optional suffix to append to end of filename
        only_continuous: bool
            whether to plot all variable types or only continuous
        '''
        print('-----Generating Plots-----')
        warnings.filterwarnings('ignore')
        self._get_cohorts()
        cohorts = ['d_f',
                   'd_m',
                   'r_f',
                   'r_m']
        if only_continuous:
            var_types = ['continuous']
        else:
            var_types = ['binary',
                         'categorical',
                         'continuous']
        for i in range(len(self.cohorts)):
            dat_to_plot = self.data[self.cohorts_bool[i]]
            for v in var_types:
                vars_to_plot = self._get_data_type(v)
                outname = v   + \
                          '_' + \
                          cohorts[i] + \
                          suffix
                filename = os.path.join(plotpath,
                                        outname)
                clarite.plot.distributions(dat_to_plot,
                                           filename=filename,
                                           continuous_kind='count',
                                           nrows=4,
                                           ncols=3,
                                           quality='low',
                                           variables=vars_to_plot)
        print('')

    def save_data(self):
        '''
        Save nhanes datasets into csv files
        '''
        print('-----Saving data-----')
        self._print_summary()
        print('')
        data_name = 'CleanNHANES_Data.csv'
        savepath  = '../Results/'
        filepath = os.path.join(savepath, data_name)
        self.data.to_csv(filepath)

        # Saving variable names
        var_lists = [self.phenotypes,
                     self.exposures]
        filenames = ['PhenotypeList.txt',
                     'ExposureList.txt']
        
        for i in range(len(var_lists)):
            output = pd.DataFrame(var_lists[i])
            output.to_csv(os.path.join(savepath, filenames[i]),
                               header=False,
                               index=False)

        # Saving weight data
        weight_data = [self.weights_discovery,
                       self.weights_replication,
                       self.weights_maps]
        filenames   = ['weights_discovery.csv',
                       'weights_replication.csv',
                       'weights_maps.csv']

        for i in range(len(weight_data)):
            weight_data[i].to_csv(os.path.join(savepath, filenames[i]))

    def _harmonize_weights_to_data(self):
        '''
        Remove weight variables that are not present in data
        and viceversa.
        Remove IDs from weight data that are not in nhanes data.

        Returns
        ----------
        data: pd.DataFrame
            data with variables without weights removed
        phenotypes: list of str
            Phenotype names without weights removed
        exposures: list of str
            Exposure names without weights removed
        weights_maps: pd.DataFrame
            Mapping data for variable - weight 
            with variables not in data removed
        weights_discovery: pd.DataFrame
            Data with weight information for the discovery
            matching the nhanes data ID
        weights_replication: pd.DataFrame
            Data with weight information for the replication
            matching the nhanes data ID
        '''
        print('-----Harmonizing weights and data-----')
        keep_phenos = list(set(self.phenotypes) & 
                           set(self.weights_maps.index))
        keep_expo   = list(set(self.exposures) & 
                           set(self.weights_maps.index))
        n_pheno_removed = len(self.phenotypes) - len(keep_phenos)
        n_expo_removed  = len(self.exposures)  - len(keep_expo)

        print('Removing ' + 
              str(n_pheno_removed) + 
              ' phenotypes, and ' + 
              str(n_expo_removed) + 
              ' exposures')

        self.phenotypes = keep_phenos
        self.exposures  = keep_expo
        self.weights_maps = self.weights_maps.loc[keep_phenos + keep_expo,]
        self._keep_data_from_lists()

        # Removing participants from nhanes data
        keep_participants = list(set(self.data.index) & 
                                 set(list(self.weights_discovery.index) +\
                                     list(self.weights_replication.index)))
        n_participants_removed = len(self.data) - len(keep_participants)

        if n_participants_removed == 0:
            print('All participants have survey weights')
        else:
            print('Removing ' + 
                  str(n_participants_removed) + 
                  ' participants')
            self.data = self.data.loc[keep_participants,]

        # Removing participants from weight data
        idx_discovery   = self.weights_discovery.index.isin(self.data.index)
        idx_replication = self.weights_replication.index.isin(self.data.index)
        self.weights_discovery   = self.weights_discovery[idx_discovery]
        self.weights_replication = self.weights_replication[idx_replication]
        print('')
        
    def _update_lists_from_data(self):
        '''
        Update phenotype and exposure lists based on data columns

        Returns
        ----------
        phenotypes: list of str
            updated phenotype list
        exposures: list of str
            updated exposure list
        '''
        for var in self.exposures + self.phenotypes:
            if var not in self.data:
                if var in self.exposures:
                    self.exposures.remove(var)
                elif var in self.phenotypes:
                    self.phenotypes.remove(var)

    def _remove_vars_in_lists(self,
                              remove_vars:list):
        '''
        Remove variables indicated by remove_vars in phenotype and
        exposure list, and call _keep_data_from_lists to update data
        
        Parameteres
        -----------
        remove_vars: list
            list of variables to remove
        
        Returns
        ----------
        phenotypes: list of str
            updated phenotype list
        exposures: list of str
            updated exposure list
        '''
        for var in remove_vars:
            if var in self.phenotypes:
                self.phenotypes.remove(var)
            elif var in self.exposures:
                self.exposures.remove(var)
        
        self._keep_data_from_lists()

    def _keep_data_from_lists(self):
        '''
        Keep data only from variables stored in covariates,
        phenotypes, and exposures lists

        Returns
        ----------
        data: pd.DataFrame
            data with only stored variables
        '''
        self.data = self.data[self.covariates + 
                              self.phenotypes + 
                              self.exposures]

    def _get_cohorts(self):
        '''
        Get cohorts based on survey year and sex

        Returns
        ----------
        cohorts: list of str
            List of cohort names
        cohorts_bool: list of bool
            List of boolean with the cohort that each participant belongs
        '''
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

    def _adjust_lipids(self):
            '''
            Adjust lipids variable if participant is on statins
            '''
            print('-----Adjusting lipid variables-----')
    
            statins = ['ATORVASTATIN_CALCIUM',
                       'SIMVASTATIN',
                       'PRAVASTATIN_SODIUM',
                       'FLUVASTATIN_SODIUM']
            for v in statins:
                if v in self.var_description:
                    print(f'\t{v} = {self.var_description[v]}')
    
            on_medication = self.data[statins].sum(axis=1) > 0
            print('There are ' +
                  str(sum(on_medication)) + 
                  ' participants on statins')
    
            self.data.loc[on_medication, 'LBDLDL'] = \
                self.data.loc[on_medication, 'LBDLDL']/0.7
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

        print('For each cohort, these are the sample sizes:')
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

    def _get_removed_vars(self,
                          data):
        '''
        Get the list of variables removed comparing data to self

        Parameters
        ----------
        data: pd.DataFrame
            data with variables removed
        
        Returns
        ----------
        removed_vars: list of str
            names of variables removed
        '''
        removed_vars = []
        for var in self.data:
            if var not in data:
                removed_vars.append(var)
        
        return(removed_vars)

    def _get_data_type(self,
                       type:str='continuous'):
        '''
        Get the list of variables that are the type

        Parameters
        ----------
        type: str
            type of variable to get

        Returns
        ----------
        var_type: list of str
            variables that are of the type
        '''
        vars = clarite.describe.get_types(self.data)
        var_type = vars[vars == type].index

        return(var_type)
