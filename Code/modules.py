import os
import pandas as pd
import numpy as np
import clarite
import warnings
import scipy.stats as stats
from functools import reduce
from statsmodels.stats.multitest import multipletests

class NhanesData:
    '''
    NHANES dataset class
    '''

    def __init__(self, dataset):
        self.data = pd.DataFrame(dataset)

    def keep_adults(self):
        '''
        Keep only adults (equal or greater than 18yo) in dataset
        '''
        print('Keeping only adults in the dataset:')
        print('-----------------------------------')
        keep_adults = self.data['RIDAGEYR'] >= 18
        self.data   = self.data.loc[keep_adults]

    def divide_discovery_replication(self):
        '''
        Separate dataset into discovery and replication based on survey cycle
        '''
        discovery   = (self.data['SDDSRVYR']==1) | (self.data['SDDSRVYR']==2)
        replication = (self.data['SDDSRVYR']==3) | (self.data['SDDSRVYR']==4)
        discovery_dat   = self.data.loc[discovery]
        replication_dat = self.data.loc[replication]
        return([discovery_dat, replication_dat])

    def divide_by_sex(self):
        '''
        Separate dataset in males and females
        '''
        cycles = self.divide_discovery_replication()

        discovery_females = NhanesData(cycles[0].loc[cycles[0]['female'] == 1])
        discovery_males   = NhanesData(cycles[0].loc[cycles[0]['male'] == 1])
        replication_females = NhanesData(cycles[1].loc[cycles[1]['female'] == 1])
        replication_males   = NhanesData(cycles[1].loc[cycles[1]['male'] == 1])
        return([discovery_females, discovery_males,
                replication_females, replication_males])

    def remove_var_without_weights(self, weights, extra_vars=[]):
        '''
        Remove variables that are not present in the survey weights
        Add extra variables that you don't want removed
        '''
        keep = set(weights.map.keys()) | set(extra_vars)
        cols = list(set(self.data.columns))

        remove_vars = []
        for i in cols:
            if i not in keep:
                print('Removing ' + i + ' variable from NHANES')
                remove_vars.append(i)

        self.data  = self.data.drop(columns=remove_vars)
        #self.data = self.data[[c for c in list(self.data) if c in keep]]

    def drop_physical_fitness(self, var_category, var_description=None):
        '''
        Drop physical fitness variables
        '''
        print('Removing physical fitness variables:')
        print('-----------------------------------')
        
        phys_fitness_vars = []
        for item in var_category.items():
            if item[1] == 'physical fitness':
                phys_fitness_vars.append(item[0])

        if var_description is not None:
            for v in phys_fitness_vars:
                print(f"\t{v} = {var_description[v]}")

        self.data = self.data.drop(columns=phys_fitness_vars)

    def drop_body_measures(self, var_category, extra_vars=[]):
        '''
        Drop body measures variables, 
        '''
        print('Removing body measures variables:')
        print('-----------------------------------')
        body_measures_vars = []
        for item in var_category.items():
            if item[1] == 'body measures' and item[0] not in extra_vars:
                body_measures_vars.append(item[0])

        self.data = self.data.drop(columns=body_measures_vars)

    def drop_disease_vars(self, var_category):
        '''
        Drop disease variables
        '''
        print('Removing disease variables:')
        print('-----------------------------------')
        disease_vars = []
        for item in var_category.items():
            if item[1] == 'disease':
                disease_vars.append(item[0])

        self.data = self.data.drop(columns=disease_vars)
        
    def drop_indeterminate_vars(self, var_description=None):
        '''
        Drop indeterminate variables
        '''
        print('Removing variables with indeterminate meanings:')
        print('-----------------------------------')
        indeterminent_vars = ["house_type","hepa","hepb","house_age",
                              "current_past_smoking","house_age","DUQ130",
                              "DMDMARTL","income","LBXSCRINV","URXUMASI",
                              "LBXSIR"]
        # Added 1/Creatinine as well
        if var_description is not None:
            for v in indeterminent_vars:
                if v in var_description:
                    print(f"\t{v} = {var_description[v]}")

        self.data = self.data.drop(columns=indeterminent_vars)

    def drop_age_groups(self):
        '''
        Drop age group categories (age is already in the covariates)
        '''
        print('Removing age categories:')
        print('-----------------------------------')
        age_groups = ["age_0_18", "age_44_65", "age_65_plus", "age_18_44"]
        self.data  = self.data.drop(columns=age_groups)

    def drop_unnecesary_vars(self, var_category):
        '''
        Call all the drop variable functions
        '''
        self.drop_physical_fitness(var_category)
        self.drop_indeterminate_vars()
        self.drop_age_groups()

    def adjust_lipids(self, var_description=None):
        '''
        Adjust lipids variable if participant is on statins
        '''
        print('Adjusting lipid variables:')
        print('-----------------------------------')

        statins = ["ATORVASTATIN_CALCIUM","SIMVASTATIN","PRAVASTATIN_SODIUM",
                    "FLUVASTATIN_SODIUM"]
        if var_description is not None:
            for v in statins:
                if v in var_description:
                    print(f"\t{v} = {var_description[v]}")

        self.data.loc[(self.data[statins].sum(axis=1) > 0), 'LBDLDL'] = \
            self.data.loc[(self.data[statins].sum(axis=1) > 0), 'LBDLDL']/0.7

    def remove_constant_vars(self, extra_vars=[]):
        '''
        Remove constant variables. Add extra_vars to include variables that you
        don't want removed (covariates, phenotypes, sex, etc)
        Clarite categorize should already been done in the dataset
        '''
        var_types    = clarite.describe.get_types(self.data)
        var_constant = var_types[var_types == 'constant'].index
        constant_variables = set(var_constant) # Get unique elements
        # Remove extra_vars from list
        for i in extra_vars:
            if i in constant_variables:
                constant_variables.remove(i)

        self.data = self.data.drop(columns=constant_variables)
        return(constant_variables)

    def remove_categorical_vars(self, extra_vars=[]):
        '''
        Remove categorical variables. Categorize should already be applied
        '''
        var_types       = clarite.describe.get_types(self.data)
        var_categorical = var_types[var_types == 'categorical'].index
        categorical_variables = set(var_categorical) # Get unique elements
        # Remove extra_vars from list
        for i in extra_vars:
            if i in categorical_variables:
                categorical_variables.remove(i)

        self.data = self.data.drop(columns=categorical_variables)
        return(categorical_variables)

    def print_unknown_vars(self, var_description):
        '''
        Print the descriptions of vars categorized as unknown
        The categorize should already be made (in remove_constant_vars)
        '''
        var_types   = clarite.describe.get_types(self.data)
        var_unknown = var_types[var_types == 'unknown'].index
        for v in list(var_unknown):
            if v in var_description:
                print(f"\t{v} = {var_description[v]}")
        return(var_unknown)

    def log_transform_skewed(self, min_max=(-0.5,0.5)):
        '''
        Log transform highly skewed continuous variables
        '''
        var_types      = clarite.describe.get_types(self.data)
        var_continuous = var_types[var_types == 'continuous'].index
        skew_test      = clarite.describe.skewness(self.data[var_continuous], dropna=True)
        skewed_vars_pos = var_continuous[skew_test['skew'] > min_max[1]]
        skewed_vars_neg = var_continuous[skew_test['skew'] < min_max[0]]
        skewed_vars     = skewed_vars_pos.append(skewed_vars_neg)

        print('Tranforming ' + str(len(skewed_vars_pos)) + 
              ' positively skewed, and ' + str(len(skewed_vars_neg)) + 
              ' negatively skewed variables')

        #Transform
        for var in skewed_vars:
            if var in skewed_vars_pos:
                vmin = self.data[var].min()
                if vmin == 0.0: # If there are zeros, add a small constant to all values
                    min_observed = np.nanmin(np.array(self.data[var])[np.nonzero(np.array(self.data[var]))])
                    c = min_observed / 10000 
                    self.data[var] = self.data[var] + c
                self.data[var] = np.log(self.data[var])
            elif var in skewed_vars_neg:
                cs = max(self.data[var]) + 1
                self.data[var] = np.log(cs - self.data[var])

    def remove_continous_outliers(self):
        '''
        Remove outliers in all continuous variables
        '''
        var_types      = clarite.describe.get_types(self.data)
        var_continuous = var_types[var_types == 'continuous'].index
        self.data = clarite.modify.remove_outliers(self.data, only=var_continuous)

    def plot_variables(self, plotpath, phenotypes, covariates, suffix=''):
        '''
        Plot variables
        '''
        os.chdir(plotpath)
    
        var_types  = clarite.describe.get_types(self.data)
        var_binary = var_types[var_types == 'binary'].index
        var_categorical = var_types[var_types == 'categorical'].index
        var_continuous  = var_types[var_types == 'continuous'].index
    
        vartypes = [var_binary, 
                    var_categorical, 
                    var_continuous, 
                    phenotypes, 
                    covariates]
        names    = ['var_binary', 
                    'var_categorical', 
                    'var_continuous', 
                    'phenotypes', 
                    'covariates']
    
        for i in range(len(vartypes)):
            outname = names[i] + suffix + '.pdf'
            clarite.plot.distributions(self.data,
                                       filename=outname,
                                       continuous_kind='count',
                                       nrows=4,
                                       ncols=3,
                                       quality='low',
                                       variables=list(vartypes[i]))

    def run_phe_ewas(self, phenotypes, covariates, survey_design):
        '''
        Run a PheEWAS across several phenotypes
        '''
        warnings.filterwarnings("ignore") # Ignore pandas warnings

        ewas_results = []
        for current_pheno in phenotypes:
            print(current_pheno)
            drop_phenos  = [p for p in phenotypes if p != current_pheno]
            data_temp = self.data.drop(columns=drop_phenos)
            data_temp = clarite.modify.rowfilter_incomplete_obs(data_temp,
                                                    only=[current_pheno])
            ewas_results_temp = clarite.analyze.ewas(outcome = current_pheno,
                                    covariates = covariates,
                                    min_n = 200,
                                    data = data_temp,
                                    regression_kind = 'weighted_glm',
                                    survey_design_spec=survey_design, 
                                    report_categorical_betas=True)
            ewas_results.append(ewas_results_temp)

        ewas_results = pd.concat(ewas_results)
        return(ewas_results)

class WeightData:
    '''
    Survey weight class
    '''

    def __init__(self, survey_weight):
        self.data = survey_weight[0] 
        self.map  = survey_weight[1]
    
    def remove_weights_not_in_data(self):
        '''
        Remove the weights variables that don't coincide between self.map and self.data
        '''
        cols = self.data.columns
        vals = set(self.map.values())

        #Which survey weights are not present
        weights_absent = []
        for i in set(vals):
            if i not in cols:
                print(i + ' weight, is not in data columns')
                weights_absent.append(i)

        #Which variables need those survey weights
        remove_vars = []
        for key, value in self.map.items():
            for i in weights_absent:
                if i == value: 
                    remove_vars.append(key)

        # Remove those vars from the weights dict
        for i in remove_vars:
            print('Removing ' + i + ' variable')
            self.map.pop(i)    

    def divide_by_nhanes(self, NhanesData1, NhanesData2):
        '''
        Divide the weight data according to the index of the NhanesData1 and 2
        '''
        survey_weight1 = WeightData([self.data, self.map])
        survey_weight2 = WeightData([self.data, self.map])
        survey_weight1.data = survey_weight1.data.loc[NhanesData1.data.index]
        survey_weight2.data = survey_weight2.data.loc[NhanesData2.data.index]
        return([survey_weight1, survey_weight2])

    def create_survey_design(self):
        '''
        Create survey desing from clarite.survey.SurveyDesignSpec
        '''
        survey_design = clarite.survey.SurveyDesignSpec(survey_df=self.data,
                                                        strata="SDMVSTRA",
                                                        cluster="SDMVPSU",
                                                        nest=True,
                                                        weights=self.map,
                                                        single_cluster='adjust',
                                                        drop_unweighted=True)
        return(survey_design)

class PheEWAS_Results:
    '''
    PheEWAS results class
    '''
    
    def __init__(self, results):
        n_results = len(results)
        # Sort indices
        for i in range(n_results):
            results[i] = results[i].sort_index()

        # Set suffixes 
        suffixes=('_df','_dm','_rf','_rm')
        for i in range(n_results):
            results[i] = results[i].add_suffix(suffixes[i])

        # Merge all datasets
        results_merged = reduce(lambda left,right: \
                         pd.merge(left, 
                                  right, 
                                  how='inner', 
                                  left_index=True, 
                                  right_index=True), 
                                  results)
        # There are some duplicated entries, so we will keep the first one
        keep = ~results_merged.duplicated()
        results_merged = results_merged[keep]

        self.results = results_merged

    def meta_analyze(self, type='total'):
        '''
        Perform a meta analysis from the pheewas results, based on the type selected
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
        
        ws = []
        for i in suffixes:
            col = 'SE' + i
            w   = np.array(1 / np.power(self.results[col], 2))
            ws.append(w)
        
        meta_se = np.sqrt( 1 / sum(ws) )

        up_term = np.zeros(meta_se.shape)
        l = 0
        for i in suffixes:
            col  = 'Beta' + i
            temp = np.array(self.results[col] * ws[l])
            up_term = up_term + temp
            l = l + 1
        
        meta_beta = up_term / sum(ws)
        zeta = meta_beta / meta_se
        pval = 2 * stats.norm.cdf(-abs(zeta))

        Ns = np.zeros(meta_se.shape)
        for i in suffixes:
            col  = 'N' + i
            temp = np.array(self.results[col])
            Ns   = Ns + temp

        self.results['SE'+out_suff]   = meta_se
        self.results['Beta'+out_suff] = meta_beta
        self.results['pvalue'+out_suff] = pval
        self.results['Variable_pvalue'+out_suff] = pval
        self.results['N'+out_suff] = Ns

    def estimate_sex_differences(self):
        '''
        Estimate sex differences, between meta analyzed male and female results
        Use suffixes to locate the female and male results
        '''
        # Estimate sex differences
        t1 = np.array(self.results['Beta_female']) - \
             np.array(self.results['Beta_male'])
        t2 = np.sqrt(np.power(np.array(self.results['SE_female']), 2) + \
                     np.power(np.array(self.results['SE_male']), 2))
        zdiff = t1 / t2
        pval  = 2*stats.norm.cdf(-abs(zdiff))
        self.results['SE_SD']   = t2
        self.results['Beta_SD'] = t1
        self.results['pvalue_SD'] = pval
        self.results['Variable_pvalue_SD'] = pval

        self.results.loc[~self.results['pvalue_SD'].isna(),
                        'pvalue_bonferroni_SD'] = multipletests(self.results.loc[~self.results['pvalue_SD'].isna(), 'pvalue_SD'],
                        method='bonferroni')[1]

    def apply_decision_tree(self):
        '''
        apply the decision tree
        '''
        #### DECISION TREE
        # 1. Z_diff is significant, also both meta female and meta male are significant, and betas are opposite
        zdiff_sign = self.results['pvalue_bonferroni_SD'] < 0.05
        both_sign  = (self.results['pvalue_female'] < 0.05 ) & \
                     (self.results['pvalue_male'] < 0.05 )
        opposite_direction = self.results['Beta_female'] * \
                             self.results['Beta_male'] < 0
        keep_qual  = zdiff_sign & \
                     both_sign & \
                     opposite_direction

        # 2. Overall nominal significance, zdiff significance bonferroni, both significant and same direction
        overall_nominal  = self.results['pvalue_total'] < 0.05
        zdiff_bonferroni = self.results['pvalue_SD'] < \
                           (0.05/sum(overall_nominal))
        same_direction   = self.results['Beta_female'] * \
                           self.results['Beta_male'] > 0
        keep_quant = overall_nominal & \
                     zdiff_bonferroni & \
                     same_direction & \
                     both_sign

        # 3. Overall nominal significance, zdiff significance boferroni, only one significant
        one_sig  = ((self.results['pvalue_female'] < 0.05 ) & \
                    (self.results['pvalue_male'] > 0.05 ) ) | \
                   ((self.results['pvalue_female'] > 0.05 ) & \
                    (self.results['pvalue_male'] < 0.05 ) )
        keep_pure = overall_nominal & \
                    zdiff_bonferroni & \
                    one_sig

        # Adding classification
        self.results['difference_type'] = 'None'
        self.results.loc[keep_qual,'difference_type']  = 'Qualitative'
        self.results.loc[keep_quant,'difference_type'] = 'Quantitative'
        self.results.loc[keep_pure,'difference_type']  = 'Pure'

    def add_variable_names(self, var_description, var_category):
        '''
        add human readable variable names, given in var_description and var_category
        '''
        index_variable  = self.results.index.get_level_values(level='Variable')
        index_phenotype = self.results.index.get_level_values(level='Outcome')
        variable_name   = _add_var_from_dict(index_variable, var_description)
        self.results['Variable_Name']  = variable_name
        phenotype_name  = _add_var_from_dict(index_phenotype, var_description)
        self.results['Outcome_Name'] = phenotype_name
        variable_category = _add_var_from_dict(index_variable, var_category)
        self.results['Variable_Category'] = variable_category
        self.results['Variable']  = index_variable
        self.results['Outcome'] = index_phenotype
        self.results = self.results.set_index(['Variable',
                                               'Variable_Name',
                                               'Outcome',
                                               'Outcome_Name'])

    def save_results(self, respath):
        '''
        Save PheEWAS results in respath
        '''
        os.chdir(respath)
        name = 'FinalResultTable.csv'
        self.results.to_csv(name)

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

def load_raw_nhanes(datapath):
    '''
    Load raw NHANES dataset
    '''
    os.chdir(os.path.join(datapath, 'nh_99-06'))
    nhanes = pd.read_csv('MainTable.csv').rename(columns={'SEQN':'ID'})\
               .set_index('ID')
    return(nhanes)

def load_clean_data(respath):
    '''
    Load cleaned data from the QC process
    '''
    os.chdir(respath)
    loadfiles = ['Discovery_Females', 'Discovery_Males', 'Replication_Females', 'Replication_Males']
    dfs = []
    for i in range(4):
        dfs.append(pd.read_csv('CleanData_' + loadfiles[i] + '.csv').set_index('ID'))
        dfs[i] = NhanesData(dfs[i])
    return(dfs)

def load_results(respath):
    '''
    Load results from the Analysis script
    '''
    os.chdir(respath)
    names = ['discovery females', 'discovery males', 'replication females', 'replication males',
             'meta females', 'meta males', 'meta total', 'sex differences']
    dfs = []
    for i in range(len(names)):
        dfs.append(pd.read_csv(names[i] + '.csv'))
    final_results = PheEWAS_Results(dfs, names)
    return(final_results)

def load_var_information(datapath, description=True):
    '''
    Load variable description file in a suitable dictionary,
    with either the variable description or the variable category
    '''
    os.chdir(os.path.join(datapath, 'nh_99-06'))
    var_description = pd.read_csv('VarDescription.csv')\
                        .drop_duplicates()\
                        .set_index('var')

    # Convert variable descriptions to a dictionary for convenience
    var_description_dict = var_description['var_desc'].to_dict()
    var_category_dict    = var_description['category'].to_dict()

    # Output required dictionary
    if description is True:
        return(var_description_dict)
    else:
        return(var_category_dict)

def load_weights(datapath, discovery=True):
    '''
    Load survey weight information,
    either from the discovery or replication dataset
    '''
    os.chdir(datapath)

    # These files map variables to their correct weights,
    # and were compiled by reading throught the NHANES codebook
    var_weights = pd.read_csv("VarWeights.csv")

    if discovery is True:
        # Read weights discovery
        survey_design_discovery = pd.read_csv("weights_discovery.txt", sep="\t")\
                                    .rename(columns={'SEQN':'ID'})\
                                    .set_index("ID")
        # Convert the data to dictionary for convenience
        weights_discovery = var_weights.set_index('variable_name')['discovery']\
                                       .to_dict()
        return([survey_design_discovery, weights_discovery])

    else:
        # Read weights replication
        survey_design_replication = pd.read_csv("weights_replication.txt", sep="\t")\
                                      .rename(columns={'SEQN':'ID'})\
                                      .set_index("ID")
        # Divide by 2 to get 4 year weights
        survey_design_replication.iloc[:,3:] = survey_design_replication.iloc[:,3:] / 2
        weights_replication = var_weights.set_index('variable_name')['replication']\
                                         .to_dict()
        return([survey_design_replication, weights_replication])

def get_phenotypes(print_descriptions=False, var_description=None, 
                   cleaned=False):
    '''
    Get the list of phenotypes and print their description if desired
    '''
    if cleaned is True:
        phenotype = ["URXUCR","LBXSCR","LBXSATSI","LBXSAL","URXUMA",
                     "LBXSAPSI","LBXSASSI","LBXSC3SI","LBXSBU","LBXBAP",
                     "LBXCPSI","LBXCRP","LBXSCLSI","LBXSCH","LBDHDL","LBDLDL",
                     "LBXSGTSI","LBXSGB","LBXGLU","LBXGH","LBXHCY",
                     "LBXSLDSI","LBXMMA","LBXSOSSI","LBXSPH","LBXSKSI",
                     "LBXSNASI","LBXSTB","LBXSCA","LBXSTP","LBXSTR","LBXSUA",
                     "LBDBANO","LBXBAPCT","LBDEONO","LBXEOPCT","LBXHCT",
                     "LBXHGB","LBDLYMNO","LBXMCHSI","LBXLYPCT","LBXMCVSI",
                     "LBXMPSI","LBDMONO","LBXMOPCT","LBXPLTSI","LBXRBCSI",
                     "LBXRDW","LBDNENO","LBXNEPCT"] # I removed the ones that were deleted in the QC process
    elif cleaned is False:
        phenotype = ["URXUCR","LBXSCR","LBXSATSI","LBXSAL","URXUMA","LBXSAPSI",
                     "LBXSASSI","LBXSC3SI","LBXSBU","LBXBAP","LBXCPSI","LBXCRP",
                     "LBXSCLSI","LBXSCH","LBDHDL","LBDLDL","LBXFER","LBXSGTSI",
                     "LBXSGB","LBXGLU","LBXGH","LBXHCY","LBXSLDSI","LBXMMA",
                     "LBXSOSSI","LBXSPH","LBXSKSI","LBXEPP","LBXSNASI","LBXTIB",
                     "LBXSTB","LBXSCA","LBXSTP","LBDPCT","LBXSTR","LBXSUA",
                     "LBDBANO","LBXBAPCT","LBDEONO","LBXEOPCT","LBXHCT",
                     "LBXHGB","LBDLYMNO","LBXMCHSI","LBXLYPCT","LBXMCVSI",
                     "LBXMPSI","LBDMONO","LBXMOPCT","LBXPLTSI","LBXRBCSI",
                     "LBXRDW","LBDNENO","LBXNEPCT","LBXIRN"]
    if print_descriptions is True and var_description is not None:
        for v in phenotype:
            if v in var_description:
                print(f"\t{v} = {var_description[v]}")

    return(phenotype)

def get_covariates():
    '''
    Get list of covariates
    '''
    covariates = ["black", "mexican", "other_hispanic", "other_eth",
                  "SES_LEVEL", "RIDAGEYR", "SDDSRVYR","BMXBMI"]
    return(covariates)

def _add_var_from_dict(index_var, var_dict):
    res = []
    for i in range(len(index_var)):
        if index_var[i] in var_dict.keys():
            res.append(var_dict[index_var[i]])
        else:
            res.append('')
    return(res)