import os
import clarite
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

from sklearn import preprocessing
from functools import reduce
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

class NhanesData:
    '''
    NHANES data class
    '''

    def __init__(self,
                 datapath:str,
                 respath:str='',
                 raw:bool=True):
        '''
        Initiates the NhanesData class and recode some values

        Parameters
        ----------
        datapath: str
            path to the raw data folder
        respath: str
            path to the results data folder 
            (only needed when not opening the raw files)
        raw: bool
            whether to load the raw or the cleaned nhanes data

        Attributes
        ----------
        data: pd.Dataframe
            Dataframe with raw nhanes dataset or cleaned data
        var_description: dict
            Dictionary with the variables description
        var_category: dict
            Dictionary with the variables category
        sex: list of str
            sex of the dataset only if loading cleaned data
        cohort: list of str
            cohort type only if loading cleaned data
        '''
        var_info  = _load_var_information(datapath)
        self.var_description = var_info[0]
        self.var_category    = var_info[1]

        if raw:
            self.data = _load_raw_nhanes(datapath)
            # SMQ077 and DDB100 have Refused/Don't Know for '7' and '9' values. 
            # We will recode them as nan
            self.data = clarite.modify.recode_values(self.data, 
                                                    {7: np.nan, 9: np.nan},
                                                    only=['SMQ077', 'DBD100'])
            #Adjust for statins (LDL=LDL/0.7, TC=TC/0.8) (2)
            self._adjust_lipids()
        else:
            self.cohort = ['Discovery', 
                           'Discovery', 
                           'Replication',
                           'Replication']
            self.sex = ['females',
                        'males',
                        'females',
                        'males']
            self._load_clean_nhanes(respath)
        
        # Divide variables
        self.divide_variables()

        print('')

    def divide_variables(self):
        '''
        Divide variables into exposures, phenotypes, and covariates

        Returns
        ----------
        dat: pd.Dataframe or list of pd.Dataframe
            dataframe with the data from covariates, phenotypes,
            and exposures
        exposures: list of str
            names of the exposure variables
        phenotypes: list of str
            names of the phenotype variables
        covariates: list of str
            names of the covariate variables
        '''
        print('-----Dividing dataset into phenotypes, exposures and covariates-----\n')
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

        # Selecting phenotypes
        phenotype_variables = []
        for i in range(len(phenotype_categories)):
            keys_list = [k for k, 
                         v in self.var_category.items()\
                         if v == phenotype_categories[i]]
            for l in keys_list:
                phenotype_variables.append(l)
        self.phenotypes = phenotype_variables

        # Selecting exposures
        exposure_variables = []
        for i in range(len(exposure_categories)):
            keys_list = [k for k, 
                         v in self.var_category.items()\
                         if v == exposure_categories[i]]
            for l in keys_list:
                exposure_variables.append(l)
        self.exposures = exposure_variables
    
        # Selecting covariates
        self.covariates = covariates

        # Generating dataset depending on whether is raw or cleaned
        if type(self.data) == pd.DataFrame:
            phenotype_dat  = self.data[phenotype_variables]
            exposure_dat   = self.data[exposure_variables]
            covariates_dat = self.data[covariates]
            # Updating data
            self.data = pd.concat([covariates_dat,
                                   phenotype_dat,
                                   exposure_dat],
                                   axis=1)
        elif type(self.data) == list:
            self._remove_not_in_data()
            for i in range(len(self.data)):
                phenotype_dat  = self.data[i][self.phenotypes]
                exposure_dat   = self.data[i][self.exposures]
                covariates_dat = self.data[i][self.covariates]
                # Updating data
                self.data[i] = pd.concat([covariates_dat,
                                          phenotype_dat,
                                          exposure_dat],
                                          axis=1)

    def keep_adults(self):
        '''
        Keep only adults (equal or greater than 18yo) in dataset

        Returns
        ----------
        data: pd.Dataframe
            data with only adults
        '''
        print('-----Keeping only adults in the dataset-----\n')
        keep_adults = self.data['RIDAGEYR'] >= 18
        print('There are only ' +
               str(sum(keep_adults)) + 
               ' adults out of ' +
               str(self.data.shape[0]) + 
               ' participants\n')

        # Keeping only adults
        self.data = self.data.loc[keep_adults]

    def remove_missing_covariates(self):
        '''
        Remove participants with missing covariates

        Returns
        ----------
        data: pd.Dataframe
            raw data with participants removed
        '''
        self.data = clarite.modify.\
                    rowfilter_incomplete_obs(self.data,
                                             only=self.covariates)

    def divide_into_cohorts(self):
        '''
        Divide the datasets into cohorts

        Returns
        ----------
        data: List of pd.Dataframe
            raw data with cohorts divided and sex removed
        cohort: List of str
            list of the cohort (either discovery or replication)
        sex: List of str
            list of the sex of the cohort (either male or female)
        covariate: List of str
            name of covariates updated without sex
        '''
        print('-----Dividing into discovery/replications and female/male cohorts-----')
        discovery   = (self.data['SDDSRVYR']==1) |\
                      (self.data['SDDSRVYR']==2) 
        replication = (self.data['SDDSRVYR']==3) |\
                      (self.data['SDDSRVYR']==4)

        female = (self.data['female'] == 1)
        male   = (self.data['male'] == 1)

        discovery_females = np.array(discovery) &\
                            np.array(female)
        replication_females = np.array(replication) &\
                              np.array(female)
        discovery_males = np.array(discovery) &\
                          np.array(male)
        replication_males = np.array(replication) &\
                            np.array(male)
        
        bools = [discovery_females,
                 discovery_males,
                 replication_females,
                 replication_males]
        new_data = []

        # Removing sex from covariates
        sexes = ['female',
                 'male']
        self.data = self.data.drop(sexes,
                                   axis=1)
        for cov in sexes:
            self.covariates.remove(cov)
        
        # Separating datasets
        for i in range(len(bools)):
            new_data.append(self.data.loc[bools[i]])

        # Setting new data
        self.data   = new_data
        self.cohort = ['Discovery', 
                       'Discovery', 
                       'Replication',
                       'Replication']
        self.sex = ['females',
                    'males',
                    'females',
                    'males']
        print('')

    def remove_missing_weights(self,
                               weights):
        '''
        Remove variables that are not present in at least one survey weight,
        except for covariates

        Returns
        ----------
        data: List of pd.Dataframe
            raw data with variables removed
        '''
        print('-----Removing variables not in the survey weights-----')
        keep = list(set(list(weights.map[0].keys()) + 
                        list(weights.map[1].keys()) +
                        self.covariates))

        cols = list(set(self.data[0].columns))

        # Removing from data
        remove_vars = []
        for i in cols:
            if i not in keep:
                remove_vars.append(i)
        for i in range(len(self.data)):
            self.data[i] = self.data[i].drop(columns=remove_vars)

        # Removing from list of variables
        self._remove_not_in_data()

    def categorize_variables(self):
        '''
        Categorize variables using clarite and remove constant variables,
        except covariates

        Returns
        ----------
        data: List of pd.Dataframe
            raw data with variables categorized and constant variables removed
        '''        
        for i in range(len(self.data)):
            print('-----Categorizing in ' + 
                  self.cohort[i] + ' ' + 
                  self.sex[i] + 
                  ' cohort-----')
            self.data[i] = clarite.modify.categorize(self.data[i])

            # Getting constant variables and removing them
            print('')
            print('-----Removing constant variables-----')
            var_types    = clarite.describe.get_types(self.data[i])
            var_constant = var_types[var_types == 'constant'].index
            constant_variables = set(var_constant)
            for cov in self.covariates:
                if cov in constant_variables:
                    constant_variables.remove(cov)
            print('We will remove ' + 
                  str(len(constant_variables)) + 
                  ' constant variables')
            self.data[i] = self.data[i].drop(columns=constant_variables)
        
        self._retain_same_across_cohorts()
        self._remove_not_in_data()

    def print_unknown_vars(self):
        '''
        Print the descriptions of vars categorized as unknown
        Remember to run categorize first
        '''
        for i in range(len(self.data)):
            print_initial = 'In the ' + \
                             self.cohort[i] + ' ' + \
                             self.sex[i] + \
                             ' cohort,'
            var_types = clarite.describe.get_types(self.data[i])
            var_unknown = var_types[var_types == 'unknown'].index
            print_end   = 'there are ' + \
                           str(len(var_unknown)) + \
                          ' variables unknown'
            print(print_initial + 
                  print_end)
            for v in list(var_unknown):
                if v in self.var_description:
                    print(f'\t{v} = {self.var_description[v]}')
                elif v not in self.var_description:
                    print(v)
            print('')

    def print_how_many_types(self):
        '''
        Print how many covariates, phenotypes, and exposures
        in the dataset are
        '''
        n_cov   = len(self.covariates)
        n_pheno = len(self.phenotypes)
        n_expo  = len(self.exposures)
        print('There are ' +
              str(n_cov) + 
              ' covariates, ' +
              str(n_pheno) + 
              ' phenotypes, and ' +
              str(n_expo) + 
              ' exposures')

    def print_how_many_participants(self):
        '''
        Print how many participants in each cohort 
        '''
        info_to_print = []
        total = 0
        for i in range(len(self.data)):
            n_participants = len(self.data[0])
            total = total + n_participants
            tp = str(n_participants) + \
                 ' in ' + \
                 str(self.cohort[i]) + ' ' + \
                 str(self.sex[i])
            info_to_print.append(tp)
        print('There are ' + 
              ', '.join(info_to_print) + 
              ' and ' + 
              str(total) + 
              ' participants')

    def make_unknown_continuous(self):
        '''
        Make unknown variables continuous

        Returns
        ----------
        data: List of pd.Dataframe
            data with unknown variables continuous
        '''
        for i in range(len(self.data)):
            var_types = clarite.describe.get_types(self.data[i])
            var_unknown = var_types[var_types == 'unknown'].index
            self.data[i] = clarite.modify.\
                           make_continuous(self.data[i],
                                           only=var_unknown)
            #self.data[i] = clarite.modify.\
            #               make_categorical(self.data[i],
            #                                only='SDDSRVYR')
        
    def make_phenotypes_continuous(self):
        '''
        Explicitly make all phenotypes continuous variables,
        which they are

        Returns
        ----------
        data: List of pd.Dataframe
            data with all phenotypes categorized as continuous
        '''
        for i in range(len(self.data)):
            self.data[i] = clarite.modify.\
                           make_continuous(self.data[i],
                                           only=self.phenotypes)

    def clarite_filters(self):
        '''
        Run clarite filters:
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
        for i in range(len(self.data)):
            print('Running in ' + 
                  self.cohort[i] + ' ' + 
                  self.sex[i] + 
                  ' cohort,')
            print('-----------------------------------')
            self.data[i] = clarite.modify.\
                           colfilter_min_n(self.data[i],
                                           skip=self.covariates)
            self.data[i] = clarite.modify.\
                           colfilter_min_cat_n(self.data[i],
                                               skip=self.covariates)
            self.data[i] = clarite.modify.\
                           colfilter_percent_zero(self.data[i],
                                                  skip=self.covariates)
        
        print('')
        self._retain_same_across_cohorts()
        self._remove_not_in_data()

    def harmonize_categorization(self):
        '''
        Compare the categorization of variables across cohorts

        Returns
        ----------
        data: List of pd.Dataframe
            data with harmonized categories
        '''
        dtypes = {}
        keys = ['one', 
                'two',
                'three',
                'four']
        for i in range(len(self.data)):
            values = clarite.describe.get_types(self.data[i])
            dtypes[keys[i]] = values

        dtypes = pd.DataFrame(dtypes)
        for i,l in dtypes.iterrows():
            m = set(l)
            m = {x for x in m if x==x}
            if len(m) > 1:
                print('There is no agreement in here ')
                print(i)
                print(l)

    def transform_continuous_log2(self):
        '''
        Log2 transform all continuous variables, and reflect those with
        negative skews

        Returns
        ----------
        data: List of pd.Dataframe
            data with log2 transformed continuous variables
        '''
        list_skewed_neg = []
        for i in range(len(self.data)):
            var_types = clarite.describe.get_types(self.data[i])
            var_continuous = var_types[var_types == 'continuous'].index
            skew_test      = clarite.describe.\
                             skewness(self.data[i][var_continuous],
                             dropna=True)
            skewed_vars_neg = var_continuous[skew_test['skew'] < 0]
            list_skewed_neg.append(skewed_vars_neg)
            print('-----Log2 transforming ' + 
                  str(len(var_continuous)) + 
                  ' variables-----')
            
        #Transform
            for var in var_continuous:
                vmin = self.data[i][var].min()
                if vmin == 0: # add a small constant to zero values
                    min_observed = np.nanmin(np.array(\
                                    self.data[i][var])[\
                                    np.nonzero(np.array(\
                                    self.data[i][var]))])
                    c = min_observed / 2
                    to_replace = self.data[i][var] == 0
                    self.data[i].loc[to_replace ,var] = c
    
            for var in var_continuous:
                #if var in skewed_vars_neg:
                #    cs = max(self.data[i][var]) + 1
                #    self.data[i][var] = np.log2(cs - 
                #                        self.data[i][var])
                #    # Reflect back
                #    cs = max(self.data[i][var]) + 0.1
                #    self.data[i][var] = cs - self.data[i][var]
                #else:
                #    self.data[i][var] = np.log2(self.data[i][var])
                self.data[i][var] = np.log2(self.data[i][var])
            
        return(list_skewed_neg)

            # Replacing columns
            #for i in range(len(self.data)):
            #    self.data[i][var_continuous] = \
            #        dat.loc[self.data[i].index, var_continuous]

    def zscore_normalization(self):
        '''
        Apply z-score normalization to values to mean center and unit variance,
        by substracting whole column by mean, and dividing by the 
        standard deviation

        Returns
        ----------
        data: List of pd.Dataframe
            data with scaled continuous variables
        '''
        for i in range(len(self.data)):
            var_types = clarite.describe.get_types(self.data[i])
            var_continuous = var_types[var_types == 'continuous'].index
            print('-----Scaling ' + 
                  str(len(var_continuous)) + 
                  ' variables-----')
            #Normalize
            self.data[i][var_continuous] = self.data[i][var_continuous].\
                                           apply(stats.zscore,
                                                 nan_policy='omit')
            #for var in var_continuous:
            #    preprocessing.scale(dat[var],
            #                        axis=0,
            #                       with_mean=True,
            #                        with_std=True,
            #                        copy=False)
            #Replace
            #for i in range(len(self.data)):
            #    self.data[i][var_continuous] = \
            #        dat.loc[self.data[i].index, var_continuous]

    def remove_continuous_outliers(self):
        '''
        Remove outliers for all continuous variables 

        Returns
        ----------
        data: list of pd.Dataframe
            data with continuous outliers removed
        '''
        for i in range(len(self.data)):
            var_types      = clarite.describe.get_types(self.data[i])
            var_continuous = var_types[var_types == 'continuous'].index
            self.data[i] = clarite.modify.\
                           remove_outliers(self.data[i],
                                           only=var_continuous)
    
    def plot_variables(self,
                       plotpath: str,
                       suffix:str = '',
                       only_continuous:bool = False):
        '''
        Plot variables

        Parameters
        ----------
        plotpath: str
            folder path to save the plots
        suffix: str
            optional suffix to append to end of
        only_continuous: bool
            whether to plot all variable types or only continuous
        '''
        warnings.filterwarnings('ignore')
        cohorts = ['d_f',
                   'd_m',
                   'r_f',
                   'r_m']
        for i in range(len(self.data)):
            var_types  = clarite.describe.get_types(self.data[i])
            var_binary = var_types[var_types == 'binary'].index
            var_categorical = var_types[var_types == 'categorical'].index
            var_continuous  = var_types[var_types == 'continuous'].index

            vartypes = [var_binary, 
                        var_categorical, 
                        var_continuous]
            names    = ['binary', 
                        'categorical', 
                        'continuous']
            if only_continuous:
                outname = names[2] + '_' + \
                          cohorts[i] + \
                          suffix
                outname = os.path.join(plotpath, outname)
                clarite.plot.distributions(self.data[i],
                                           filename=outname,
                                           continuous_kind='count',
                                           nrows=4,
                                           ncols=3,
                                           quality='low',
                                           variables=list(vartypes[2]))
            elif only_continuous == False:
                for l in range(len(vartypes)):
                    if len(vartypes[l]) > 0:
                        outname = names[l] + '_' + \
                                  cohorts[i] + \
                                  suffix
                        outname = os.path.join(plotpath, outname)
                        clarite.plot.distributions(self.data[i],
                                                   filename=outname,
                                                   continuous_kind='count',
                                                   nrows=4,
                                                   ncols=3,
                                                   quality='low',
                                                   variables=list(vartypes[l]))

    def save_data(self,
                  savepath:str):
        '''
        Save nhanes datasets into csv files

        Parameters
        ----------
        savepath: str
            path to save the csv files
        '''
        # Order columns
        self._order_columns()
        for i in range(len(self.data)):
            name = 'CleanData_' + \
                   self.cohort[i] + '_' + \
                   self.sex[i] + \
                   '.csv'
            filepath = os.path.join(savepath, name)
            self.data[i].to_csv(filepath)

    def save_pheno_expo_list(self,
                             savepath:str):
        '''
        Save file with list of phenotypes and exposures

        Parameters
        ----------
        savepath: str
            path to save the csv files
        '''
        filepaths = [os.path.join(savepath, 'Phenotypes.txt'),
                     os.path.join(savepath, 'Exposures.txt')]
        name_lists = [self.phenotypes,
                      self.exposures]
        
        for i in range(len(filepaths)):
            variable_description = _add_var_from_dict(name_lists[i],
                                                      self.var_description)
            output = pd.DataFrame(name_lists[i])
            output['Description'] = variable_description
            output.to_csv(filepaths[i],
                          header=False,
                          index=False)

    def run_phe_ewas(self,
                     weights):
        '''
        Run a PheEWAS

        Parameters
        ----------
        weights: WeightData
            weights harmonized and with survey designs created

        Returns
        ----------
        results: PheEWAS_Results
            results from PheEWAS
        '''
        warnings.filterwarnings('ignore') # Ignore pandas warnings

        total_results = []
        for i in range(len(self.data)):
            res_temp = clarite.analyze.association_study(data=self.data[i],
                                           outcomes=self.phenotypes,
                                           covariates=self.covariates,
                                           regression_kind='weighted_glm',
                                           survey_design_spec=weights.design[i],
                                           report_categorical_betas=True)
            total_results.append(res_temp)
        results = PheEWAS_Results(total_results)
        return(results)

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

    def _remove_variables_from_list(self,
                                    to_remove:list,
                                    to_print:bool):
        '''
        Remove variables from list of exposures and phenotypes

        Parameters
        ----------
        to_remove: list
            name of variables to remove
        print: bool
            whether to print the variables and types that are
            being removed
        '''
        for var in to_remove:
            if var in self.exposures:
                self.exposures.remove(var)
                if to_print:
                    print('Removing ' + 
                          var + 
                          ' from exposures')
            elif var in self.phenotypes:
                self.phenotypes.remove(var)
                if to_print:
                    print('Removing ' + 
                          var + 
                          ' from phenotypes')
            else:
                print(var + 
                      ' not in phenotypes or exposures')
        print('')

    def _remove_not_in_data(self):
        '''
        Remove variable names in exposures and phenotypes
        not seen in the dataframe

        Returns
        ----------
        exposures: list of str
            name of exposures matching data
        phenotypes: list of str
            name of phenotypes matching data
        covariates: list of str
            name of covariates matching data
        '''
        print('-----Matching variable names with data-----\n')
        self.phenotypes = _remove_setdiff(self.phenotypes, 
                                          self.data[0].columns, 
                                          'phenotype')
        self.exposures = _remove_setdiff(self.exposures, 
                                         self.data[0].columns,
                                         'exposure')
        self.covariates = _remove_setdiff(self.covariates,
                                          self.data[0].columns,
                                          'covariate')
        print('')

    def _retain_same_across_cohorts(self):
        '''
        Returns list of variables that are common across cohorts

        Returns
        ----------
        data: List of pd.Dataframe
            list of dataframes with the same variables
        '''
        print('-----Keeping the same variables across cohorts-----')
        for i in range(len(self.data)):
            s = set(list(self.data[i]))
            if i == 0:
                all_vars = s
            else:
                all_vars = all_vars & s

        for i in range(len(self.data)):
            self.data[i] = self.data[i][all_vars]

        print(f'{len(all_vars)} variables in common\n')

    def _order_columns(self):
        '''
        Order the datasets columns by covariates, phenotypes, and exposures

        Returns
        ----------
        data: list of pd.Dataframe
            datasets with columns ordered
        '''
        new_data = []
        for i in range(len(self.data)):
            dat = pd.concat([self.data[i][self.covariates],
                            self.data[i][self.phenotypes],
                            self.data[i][self.exposures]],
                            axis=1)
            new_data.append(dat)
        
        self.data = new_data

    def _load_clean_nhanes(self,
                           respath:str):
        '''
        Load cleaned nhanes data resulted from the QC process

        Parameters
        ----------
        respath: str
            folder path to the cleaned data files
        '''
        dfs = []
        for i in range(len(self.cohort)):
            filename = 'CleanData_' + \
                       self.cohort[i] + '_' +\
                       self.sex[i] + \
                       '.csv'
            filepath = os.path.join(respath,
                                    filename)
            dfs.append(pd.read_csv(filepath).set_index('ID'))
        self.data = dfs

class WeightData:
    '''
    Survey weight class
    '''

    def __init__(self,
                 datapath):
        '''
        Initiates the WeightData class

        Parameters
        ----------
        datapath: str
            Str with the path to the data folder

        Attributes
        ----------
        data: List of Dataframe
            List of data with weight values per participant
        map: List of dict
            List of dict with the maping of variables and weight variable
        cohort: List of str
            List of cohort ('Discovery' or 'Replication')
        '''
        weights_discovery   = _load_weights(datapath)
        weights_replication = _load_weights(datapath, False)

        self.data = [weights_discovery[0], weights_replication[0]]
        self.map  = [weights_discovery[1], weights_replication[1]]
        self.cohort = ['Discovery', 'Replication']
    
    def harmonize(self):
        '''
        Remove the weights variables that don't coincide between map and data
        '''
        print('-----Harmonizing survey weights-----')
        for c in range(len(self.data)):
            cols = self.data[c].columns
            vals = set(self.map[c].values())

            # Which survey weights are not present
            weights_absent = []
            for i in set(vals):
                if i not in cols:
                    print(i + 
                          ' weight, is not in data columns in ' +
                          self.cohort[c] + 
                          ' cohort')
                    weights_absent.append(i)

            # Which variables need those survey weights
            remove_vars = []
            for key, value in self.map[c].items():
                for i in weights_absent:
                    if i == value:
                        remove_vars.append(key)

            # Remove those vars from the weights dict
            for i in remove_vars:
                print('Removing ' + 
                      i + 
                      ' variable in ' +
                      self.cohort[c] + 
                      ' cohort')
                self.map[c].pop(i)
        print('')

    def harmonize_with_nhanes(self, 
                              nhanes: NhanesData):
        '''
        Keep only the IDs that are present in the nhanes data and
        split in the same way as nhanes

        Parameters
        ----------
        nhanes: NhanesData
            object with the nhanes data divided into cohorts
        
        Returns
        ----------
        data: pd.Dataframe
            data matching the IDs from nhanes 
        '''
        print('-----Harmonizing weights with Nhanes data-----')
        new_data = []
        for i in range(len(nhanes.cohort)):
            indice = [a for a, 
                      x in enumerate(self.cohort)\
                      if x==nhanes.cohort[i]]
            n = self.data[indice[0]].loc[\
                    nhanes.data[i].index]
            new_data.append(n)
        self.data = new_data
        print('')
        
    def create_survey_design(self):
        '''
        Create survey desing from clarite.survey.SurveyDesignSpec.
        Use after harmonizing with nhanes

        Returns
        ----------
        design = list of SurveyDesignSpec
            survey design
        '''
        print('-----Creating survey designs-----')
        survey_design = []
        for i in range(len(self.data)):
            if i < 2:
                indice = 0
            else:
                indice = 1
            s = clarite.survey.SurveyDesignSpec(survey_df=self.data[i],
                                                strata='SDMVSTRA',
                                                cluster='SDMVPSU',
                                                nest=True,
                                                weights=self.map[indice],
                                                single_cluster='adjust',
                                                drop_unweighted=True)
            survey_design.append(s)
        self.design = survey_design
        print('')

class PheEWAS_Results:
    '''
    PheEWAS results class
    '''
    
    def __init__(self, 
                 results:list):
        '''
        Initialize class

        Attributes
        ----------
        data: pd.Dataframe
            result data from PheEWAS analysis
        suffixes: list of str
            names to distinguish between different results
        N: list of int
            number of tests in each category
        converged: list of int
            number of tests that converged
        not_converged: list of int
            number of tests that didn't converged
        '''
        n_results = len(results)
        # Sort indices
        for i in range(n_results):
            results[i] = results[i].sort_index()

        # Set N and (not)converged
        self.N = []
        self.converged=[]
        self.not_converged=[]
        for i in range(n_results):
            nrows = results[i].shape[0]
            nc  = sum(results[i]['Converged'])
            nnc = sum(results[i]['Converged']==False)
            self.N.append(nrows)
            self.converged.append(nc)
            self.not_converged.append(nnc)

        # Set suffixes 
        self.suffixes=['_df','_dm','_rf','_rm']
        for i in range(n_results):
            results[i] = results[i].add_suffix(self.suffixes[i])
        
        # Merge all datasets
        results_merged = reduce(lambda left,right: \
                         pd.merge(left, 
                                  right, 
                                  how='inner', 
                                  left_index=True, 
                                  right_index=True), 
                                  results)
        # There are some duplicated entries, with the same values
        # so we will keep the first one
        keep = ~results_merged.duplicated()
        results_merged = results_merged[keep]

        self.data = results_merged

    def meta_analyze(self, 
                     type:str='total'):
        '''
        Perform a meta analysis from the PheEWAS results, 
        based on the type selected

        Parameters
        ----------
        type: str
            Type of meta analysis to generate ('male', 'female', or 'total')

        Returns
        ----------
        data: pd.Dataframe
            updated dataframe with meta analysis values
        suffixes: list of str
            updated list of suffixes
        '''
        if type == 'total':
            suffixes = ('_female','_male')
            out_suff = ('_total')
        elif type == 'female':
            suffixes = ('_df', '_rf')
            out_suff = ('_female')
        elif type == 'male':
            suffixes = ('_dm', '_rm')
            out_suff = ('_male')
        
        self.suffixes.append(out_suff)
        ws = []
        for i in suffixes:
            col = 'SE' + i
            w   = np.array(1 / np.power(self.data[col], 2))
            ws.append(w)
        
        meta_se = np.sqrt( 1 / sum(ws) )

        up_term = np.zeros(meta_se.shape)
        l = 0
        for i in suffixes:
            col  = 'Beta' + i
            temp = np.array(self.data[col] * ws[l])
            up_term = up_term + temp
            l = l + 1
        
        meta_beta = up_term / sum(ws)
        zeta = meta_beta / meta_se
        pval = 2 * stats.norm.cdf(-abs(zeta))

        Ns = np.zeros(meta_se.shape)
        for i in suffixes:
            col  = 'N' + i
            temp = np.array(self.data[col])
            Ns   = Ns + temp

        self.data['SE'+out_suff]   = meta_se
        self.data['Beta'+out_suff] = meta_beta
        self.data['pvalue'+out_suff] = pval
        self.data['Variable_pvalue'+out_suff] = pval
        self.data['N'+out_suff] = Ns

    def estimate_sex_differences(self):
        '''
        Estimate sex differences, between meta analyzed male and female results
        Use suffixes to locate the female and male results

        Returns
        ----------
        data: pd.Dataframe
            updated results dataframe
        '''
        # Estimate sex differences
        t1 = np.array(self.data['Beta_female']) - \
             np.array(self.data['Beta_male'])
        t2 = np.sqrt(np.power(np.array(self.data['SE_female']), 2) + \
                     np.power(np.array(self.data['SE_male']), 2))
        zdiff = t1 / t2
        pval  = 2*stats.norm.cdf(-abs(zdiff))
        self.data['SE_SD']   = t2
        self.data['Beta_SD'] = t1
        self.data['pvalue_SD'] = pval
        self.data['Variable_pvalue_SD'] = pval

        self.data.loc[~self.data['pvalue_SD'].isna(),
                      'pvalue_bonferroni_SD'] = multipletests(self.data.loc[~self.data['pvalue_SD'].isna(), 'pvalue_SD'],
                      method='bonferroni')[1]

    def apply_decision_tree(self):
        '''
        Apply the decision tree base on Winkler et al 2017 without any
        previous hypothesis.
        It will select all Bonferroni corrected associations that are
        different between sexes, and apply a filtering by overall association
        as well.

        Returns
        ----------
        data: pd.Dataframe
            updated dataframe with the type of difference
        '''
        print('-----Applying decision criteria-----')
        #### Total difference
        total_diff = self.data['pvalue_bonferroni_SD'] < 0.05

        #### Filtering first
        filter_t   = 10 ** -5
        overall_filter = self.data['pvalue_total'] < filter_t
        bonf_t = 0.05 / sum(overall_filter)
        filter_diff = self.data['pvalue_SD'] < bonf_t

        significants = total_diff | \
                       filter_diff    

        #### CLASSIFICATION
        # 1. Sinificant, both meta female and meta male are significant, and betas are opposite
        both_sign  = (self.data['pvalue_female'] < 0.05 ) & \
                     (self.data['pvalue_male'] < 0.05 )
        opposite_direction = self.data['Beta_female'] * \
                             self.data['Beta_male'] < 0
        keep_qual  = significants & \
                     both_sign & \
                     opposite_direction

        # 2. Overall nominal significance, zdiff significance bonferroni, both significant and same direction
        same_direction   = self.data['Beta_female'] * \
                           self.data['Beta_male'] > 0
        keep_quant = significants & \
                     same_direction & \
                     both_sign

        # 3. Overall nominal significance, zdiff significance boferroni, only one significant
        one_sig  = ((self.data['pvalue_female'] < 0.05 ) & \
                    (self.data['pvalue_male'] > 0.05 ) ) | \
                   ((self.data['pvalue_female'] > 0.05 ) & \
                    (self.data['pvalue_male'] < 0.05 ) )
        keep_pure = significants & \
                    one_sig

        print('There are ' + 
              str(sum(significants)) +
              ' significant results, ' + 
              str(sum(keep_qual)) + 
              ' qualitative, ' +
              str(sum(keep_quant)) +
              ' quantitative, and ' +
              str(sum(keep_pure)) +
              ' pure')    

        # Adding classification
        self.data['difference_type'] = 'None'
        self.data.loc[keep_qual,'difference_type']  = 'Qualitative'
        self.data.loc[keep_quant,'difference_type'] = 'Quantitative'
        self.data.loc[keep_pure,'difference_type']  = 'Pure'

    def add_variable_names(self, 
                           var_description:dict, 
                           var_category:dict):
        '''
        add human readable variable names, 
        given in var_description and var_category
        
        Parameters
        ----------
        var_description: dict
            description of nhanes variable name
        var_category: dict
            category of nhanes variable name

        Returns
        ----------
        data: pd.Dataframe
            dataframe with nhanes variable names expanded
        '''
        index_variable  = self.data.index.get_level_values(level='Variable')
        index_phenotype = self.data.index.get_level_values(level='Outcome')
        variable_name   = _add_var_from_dict(index_variable, var_description)
        self.data['Variable_Name']  = variable_name
        phenotype_name  = _add_var_from_dict(index_phenotype, var_description)
        self.data['Outcome_Name'] = phenotype_name
        variable_category = _add_var_from_dict(index_variable, var_category)
        self.data['Variable_Category'] = variable_category
        self.data['Variable']  = index_variable
        self.data['Outcome'] = index_phenotype
        self.data = self.data.set_index(['Variable',
                                         'Variable_Name',
                                         'Outcome',
                                         'Outcome_Name'])

    def plot_partial_regression(self,
                                dat_f,
                                dat_m,
                                covariates:list):
        '''
        Plot scatter plots from sex different significant 
        Variable-Outcome regressions.
        You need to have ran estimate_sex_differences(), 
        and apply_decision_tree()

        Parameters
        ----------
        dat_f: pd.DataFrame
            dataframe with the female cohort
        dat_m: pd.DataFrame
            dataframe with the male cohort
        covariates: list of str
            name of the covariates to control for
        '''
        not_none = self.data['difference_type'] != 'None'
        dat = pd.concat([dat_f,
                         dat_m])
        for row in self.data[not_none].iterrows():
            predictors = covariates + [row[0][1]]
            dat_temp   = dat_f[predictors + [row[0][0]]]
            dat_temp   = dat_temp.dropna()
            var_residuals_f = _get_residuals(dat_temp,
                                             predictors,
                                             row[0][0])
            predictors = covariates + [row[0][0]]
            outcome_residuals_f = _get_residuals(dat_temp,
                                                 predictors,
                                                 row[0][1])

        ##IM DOING THIS!!!!!!!!!!!

    def save_results(self, 
                     respath:str):
        '''
        Save PheEWAS results in respath

        Parameters
        ----------
        respath: str
            folder path to save results file
        '''
        name = 'FinalResultTable.csv'
        filepath = os.path.join(respath, name)
        self.data.to_csv(filepath)

def set_project_paths():
    '''
    Setting paths for the project
    Outputs the main, data and results paths
    '''
    os.chdir('..')
    mainpath = os.getcwd()
    datapath = os.path.join(mainpath, 'Data')
    respath  = os.path.join(mainpath, 'Results')
    return([mainpath, datapath, respath])

def _remove_setdiff(list1:list,
                    list2:list,
                    type:str):
    '''
    Remove elements in list1 not in list2

    Parameters
    ----------
    list1: list
    list2: list
    type: str
        type of list to print

    Returns
    ----------
    list1: list
        list1 with elements not seen in list2 removed
    '''
    remove = np.setdiff1d(list1,
                          list2)
    for el in remove:
        print('Removing ' +
              el + ' ' +
              type)
        list1.remove(el)
    return(list1)

def _load_raw_nhanes(datapath:str):
    '''
    Load raw NHANES dataset

    Parameters
    ----------
    datapath: str
       path with the folder to the dataset
    '''
    filename = os.path.join(datapath, 'nh_99-06', 'MainTable.csv')
    nhanes   = pd.read_csv(filename).\
                  rename(columns={'SEQN':'ID'})\
                  .set_index('ID')
    return(nhanes)

def _load_var_information(datapath:str):
    '''
    Load variable description file in a suitable dictionary,
    with either the variable description or the variable category

    Parameters
    ----------
    datapath: str
        folder path for the data location
    '''
    filename = os.path.join(datapath, 'nh_99-06', 'VarDescription.csv')
    var_description = pd.read_csv(filename)\
                        .drop_duplicates()\
                        .set_index('var')

    # Convert variable descriptions to a dictionary for convenience
    var_description_dict = var_description['var_desc'].to_dict()
    var_category_dict    = var_description['category'].to_dict()

    # Output list of dict
    return([var_description_dict,
            var_category_dict])

def _load_weights(datapath, discovery=True):
    '''
    Load survey weight information,
    either from the discovery or replication dataset
    '''
    filename_var = os.path.join(datapath,
                                'VarWeights.csv')
    filename_disc = os.path.join(datapath,
                                 'weights_discovery.txt')
    filename_repl = os.path.join(datapath,
                                 'weights_replication.txt')

    # These files map variables to their correct weights,
    # and were compiled by reading throught the NHANES codebook
    var_weights = pd.read_csv(filename_var)

    if discovery is True:
        # Read weights discovery
        survey_design_discovery = pd.read_csv(filename_disc, 
                                              sep='\t')\
                                    .rename(columns={'SEQN':'ID'})\
                                    .set_index('ID')
        # Convert the data to dictionary for convenience
        weights_discovery = \
            var_weights.set_index('variable_name')\
                ['discovery'].to_dict()

        return([survey_design_discovery,
                weights_discovery])

    else:
        # Read weights replication
        survey_design_replication = pd.read_csv(filename_repl,
                                                sep='\t')\
                                      .rename(columns={'SEQN':'ID'})\
                                      .set_index('ID')
        # Divide by 2 to get 4 year weights
        survey_design_replication.iloc[:,3:] = \
            survey_design_replication.iloc[:,3:] / 2

        weights_replication = \
            var_weights.set_index('variable_name')\
                ['replication'].to_dict()

        return([survey_design_replication,
                weights_replication])

def _add_var_from_dict(index_var, var_dict):
    res = []
    for i in range(len(index_var)):
        if index_var[i] in var_dict.keys():
            res.append(var_dict[index_var[i]])
        else:
            res.append('')
    return(res)

def _get_residuals(dat,
                   predictors:list,
                   outcome:str):
    '''
    Get the residuals from a linear regression

    Parameters
    ----------
    dat: pd.DataFrame
        data with no NAs
    predictors: list of str
        name of the predictors for the regression
    outcome: str
        name of the outcome
    '''
    regression = LinearRegression().\
                  fit(dat[predictors], dat[outcome])
    residuals  = dat[outcome] - \
                  regression.\
                  predict(dat[predictors])
    return(residuals)

def check_balanced_tests(dat,
                         savepath:str):
    '''
    For every exposure-phenotype test, check the precentage of females
    and output the result

    Parameters
    ----------
    dat: pd.DataFrame
        results from the PheEWAS analysis
    savepath: str
        path to save the list of percentages
    '''
    size_names = ['N_df',
                  'N_dm',
                  'N_rf',
                  'N_rm']
    
    N_f = dat[size_names[0]] + dat[size_names[2]]
    N_m = dat[size_names[1]] + dat[size_names[3]]
    N_t = N_f + N_m
    percent = N_f / N_t * 100
    percent.round(2)
    new_dat = pd.DataFrame(percent.round(2), 
                           columns=['percent'])
    name = 'PercentFemales.csv'
    filepath = os.path.join(savepath,
                            name)
    new_dat.to_csv(filepath)
    