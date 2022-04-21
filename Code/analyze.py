import os
import utils
import random
import weights
import clarite
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats

from functools import reduce
from statsmodels.stats.multitest import multipletests

class NhanesClean:
    '''
    Nhanes cleaned class
    '''

    def __init__(self):
        '''
        Initiate the class

        Attributes
        ----------
        data: pd.DataFrame
            NHANES clean data
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
        print('-----Loading NHANES cleaned data-----')

        # Setting data path
        data_path = '../Results/'

        # Loading data
        var_info  = utils.load_var_information()
        self.var_description = var_info[0]
        self.var_category    = var_info[1]

        filename = os.path.join(data_path,
                                'CleanNHANES_Data.csv')
        nhanes   = pd.read_csv(filename).\
                      set_index('ID')
        self.data = nhanes

        # Dividing data
        self.covariates = utils.get_covariate_list()
        expos = pd.read_csv(os.path.join(data_path,
                                         'ExposureList.txt'),
                            header=None)
        self.exposures = expos[0].values.tolist()

        phenos = pd.read_csv(os.path.join(data_path,
                                          'PhenotypeList.txt'),
                             header=None)
        self.phenotypes = phenos[0].values.tolist()

        # Weight data
        weight = weights.load_weights_cleaned()
        self.weights_discovery   = weight[0]
        self.weights_replication = weight[1]
        self.weights_maps        = weight[2]

        # Cohorts
        self._get_cohorts()

    def categorize_variables(self):
        '''
        Categorize variables using clarite

        Returns
        ----------
        data: pd.DataFrame
            raw data with variables categorized
        '''
        print('-----Categorizing data-----')
        self.data = clarite.modify.categorize(self.data)

        var_types    = clarite.describe.get_types(self.data)
        var_unknown  = var_types[var_types == 'unknown'].index
        if len(var_unknown) > 0:
            self.data = clarite.modify.\
                                make_continuous(self.data,
                                                only=var_unknown)
        print('') 

    def create_survey_design(self):
        '''
        Create survey desing from clarite.survey.SurveyDesignSpec

        Returns
        ----------
        weights_design = list of SurveyDesignSpec
            survey designs following the cohorts division
        '''
        print('-----Creating survey designs-----')
        survey_design = []
        for i in range(len(self.cohorts)):
            keep_idx = self.data[self.cohorts_bool[i]].index
            if 'discovery' in self.cohorts[i]:
                dat = self.weights_discovery.loc[keep_idx,:]
                weight_map = self.weights_maps['discovery'].to_dict()
            elif 'replication' in self.cohorts[i]:
                dat = self.weights_replication.loc[keep_idx,:]
                weight_map = self.weights_maps['replication'].to_dict()

            s = clarite.survey.SurveyDesignSpec(survey_df=dat,
                                            strata='SDMVSTRA',
                                            cluster='SDMVPSU',
                                            nest=True,
                                            weights=weight_map,
                                            single_cluster='adjust',
                                            drop_unweighted=True)
            survey_design.append(s)
        self.weights_design = survey_design
        print('')

    def run_phe_ewas(self):
        '''
        Run a PheEWAS evaluating all phenotypes with all exposures
        while controlling for covariates

        Returns
        ----------
        results: PheEWAS_Results
            results from PheEWAS
        '''
        warnings.filterwarnings('ignore') # Ignore pandas warnings

        total_results = []
        covs = self.covariates.copy()
        remove_sex = ['female',
                      'male']
        for c in remove_sex:
            covs.remove(c)

        for i in range(len(self.cohorts)):
            dat  = self.data[self.cohorts_bool[i]]
            dat  = dat[covs + 
                       self.phenotypes +
                       self.exposures]
            res_temp = clarite.analyze.association_study(data=dat,
                                    outcomes=self.phenotypes,
                                    covariates=covs,
                                    regression_kind='weighted_glm',
                                    survey_design_spec=self.weights_design[i],
                                    report_categorical_betas=True)
            total_results.append(res_temp)
        results = PHE_EWAS(total_results)
        return(results)

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

    def _reduce_variables(self,
                          n_phenotypes:int=10,
                          n_exposures:int=10):
        '''
        Reduce the number of exposures and phenotypes to
        a random set.
        Used for testing analyses and not wait for the 
        entire dataset

        Parameters
        ----------
        n_phenotypes: int
            number of phenotypes to retain
        n_exposures: int
            number of exposures to retain
        '''

        self.exposures = random.sample(self.exposures,
                                       n_exposures)
        self.phenotypes = random.sample(self.phenotypes,
                                        n_phenotypes)
        
        self.data = self.data[self.covariates + 
                              self.phenotypes + 
                              self.exposures]

class PHE_EWAS:
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
        converged: list of int
            number of tests that converged
        not_converged: list of int
            number of tests that didn't converged
        '''
        n_results = len(results)
        # Sort indices
        for i in range(n_results):
            results[i] = results[i].sort_index()

        # Set (not)converged
        self.converged=[]
        self.not_converged=[]
        for i in range(n_results):
            nrows = results[i].shape[0]
            nc  = sum(results[i]['Converged'])
            nnc = sum(results[i]['Converged']==False)
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
        
        # Keep only converged results
        converged_cols = ['Converged_df',
                          'Converged_rf',
                          'Converged_dm',
                          'Converged_rm']
        keep = results_merged[converged_cols].all(axis=1)

        self.data = results_merged[keep]

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

        self.data['SE' + out_suff]   = meta_se
        self.data['Beta' + out_suff] = meta_beta
        self.data['pvalue' + out_suff] = pval
        self.data['Variable_pvalue' + out_suff] = pval
        self.data['N' + out_suff] = Ns

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
        # 1. Significant, both meta female and meta male are significant, and betas are opposite
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
        variable_name   = utils.add_var_from_dict(index_variable,
                                                   var_description)
        self.data['Variable_Name']  = variable_name
        phenotype_name  = utils.add_var_from_dict(index_phenotype,
                                                  var_description)
        self.data['Outcome_Name'] = phenotype_name
        variable_category = utils.add_var_from_dict(index_variable,
                                                    var_category)
        self.data['Variable_Category'] = variable_category
        self.data['Variable']  = index_variable
        self.data['Outcome'] = index_phenotype
        self.data = self.data.set_index(['Variable',
                                         'Variable_Name',
                                         'Outcome',
                                         'Outcome_Name'])

    def save_results(self,
                     respath: str = '../Results/'):
        '''
        Save PheEWAS results in respath

        Parameters
        ----------
        respath: str
            folder path to save results file
        '''
        name = 'CompleteResultsTable.csv'
        filepath = os.path.join(respath, name)
        self.data.to_csv(filepath)

        # Saving significant only
        name = 'SignificantResultsTable.csv'
        sig_bool = self.data.loc[:, 'difference_type'] != 'None'
        sig_table = self.data.loc[sig_bool]
        filepath = os.path.join(respath, name)
        sig_table.to_csv(filepath)

        # Creating tables
        columns_keep = ['Variable_Category',
                        'Variable_Name',
                        'Outcome_Name',
                        'Beta_female',
                        'pvalue_female',
                        'Beta_male',
                        'pvalue_male',
                        'pvalue_SD']
        columns_round = ['Beta_female',
                         'Beta_male']
        formats = {'pvalue_female': '{:.2E}',
                   'pvalue_male': '{:.2E}',
                   'pvalue_SD': '{:.2E}'}
        types = ['Pure',
                 'Quantitative',
                 'Qualitative']
        data = self.data.reset_index()
        for t in types:
            name = t + 'Table.csv'
            filepath = os.path.join(respath, name)
            table_bool = data.loc[:, 'difference_type'] == t
            table = data.loc[table_bool, columns_keep]
            for col, f in formats.items():
                table[col] = table[col].map(lambda x: f.format(x))
            table[columns_round] = table[columns_round].round(3)
            # Sort values
            table = table.sort_values(by=['Variable_Category',
                                          'Variable_Name',
                                          'Outcome_Name'])
            table.to_csv(filepath, sep='&', index=False)
