import os
import pandas as pd
import numpy as np
import clarite
import warnings
import scipy.stats as stats
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
        keep = set(weights[1].keys()) | set(extra_vars)

        self.data = self.data[[c for c in list(self.data) if c in keep]]

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
                              "DMDMARTL","income"]
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

        print('Tranforming ' + str(len(skewed_vars_pos)) + ' positevely skewed, and ' +
               str(len(skewed_vars_neg)) + ' negatevely skewed variables')

        #Transform
        for var in skewed_vars:
            if var in skewed_vars_pos:
                vmin = self.data[var].min()
                if vmin == 0.0: # If there are zeros, add a small constant to all values
                    min_observed = np.nanmin(np.array(self.data[var])[np.nonzero(np.array(self.data[var]))])
                    c = min_observed / 100 
                    self.data[var] = self.data[var] + c
                self.data[var] = np.log(self.data[var])
            elif var in skewed_vars_neg:
                cs = max(self.data[var]) + 1
                self.data[var] = np.log(cs - self.data[var])

    def plot_variables(self, plotpath, phenotypes, covariates, suffix=''):
        '''
        Plot variables
        '''
        os.chdir(plotpath)
    
        var_types  = clarite.describe.get_types(self.data)
        var_binary = var_types[var_types == 'binary'].index
        var_categorical = var_types[var_types == 'categorical'].index
        var_continuous  = var_types[var_types == 'continuous'].index
    
        vartypes = [var_binary, var_categorical, var_continuous, phenotypes, covariates]
        names    = ['var_binary', 'var_categorical', 'var_continuous', 'phenotypes', 'covariates']
    
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
            ewas_results_temp = clarite.analyze.ewas(phenotype = current_pheno,
                                                     covariates = covariates,
                                                     min_n = 200,
                                                     data = self.data.drop(columns=drop_phenos),
                                                     regression_kind = 'weighted_glm',
                                                     survey_design_spec=survey_design)
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
                print(i + ' is not in cols')
                weights_absent.append(i)

        #Which variables need those survey weights
        remove_vars = []
        for key, value in self.map.items():
            for i in weights_absent:
                if i == value: 
                    remove_vars.append(key)

        # Remove those vars from the weights dict
        for i in remove_vars:
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
                                                        single_cluster='adjust')
        return(survey_design)

class PheEWAS_Results:
    '''
    PheEWAS results class
    '''
    
    def __init__(self, results, name_of_results):
        n_results = len(results)
        # Sort based on first dataframe
        for i in range(1, n_results):
            results[i] = results[i].reindex(results[0].index)
        self.results = results
        self.name    = name_of_results

    def meta_analyze(self, name, indices=(0,1,2,3)):
        '''
        Perform a meta analysis from the pheewas results, based on the indices selected
        '''
        ws = []
        for i in indices:
            w = np.array(1 / np.power(self.results[i]['SE'], 2))
            ws.append(w)
        
        meta_se = np.sqrt( 1 / sum(ws) )

        up_term = np.zeros(meta_se.shape)
        l = 0
        for i in indices:
            temp = np.array(self.results[i]['Beta'] * ws[l])
            up_term = up_term + temp
            l = l + 1
        
        meta_beta = up_term / sum(ws)
        zeta = meta_beta / meta_se
        pval = 2 * stats.norm.cdf(-abs(zeta))

        Ns = np.zeros(meta_se.shape)
        for i in indices:
            temp = np.array(self.results[i]['N'])
            Ns   = Ns + temp

        mindex   = self.results[0].index
        final_df = pd.DataFrame(index=mindex)
        final_df['SE']   = meta_se
        final_df['Beta'] = meta_beta
        final_df['pvalue'] = pval
        final_df['Variable_pvalue'] = pval
        final_df['N'] = Ns

        self.results.append(final_df)
        self.name.append(name)

    def estimate_sex_differences(self, indices=(4,5)):
        '''
        Estimate sex differences, between meta analyzed male and female results
        Use indices to locate the female and male results
        '''
        # Estimate sex differences
        females_i = indices[0]
        males_i   = indices[1]
        t1 = np.array(self.results[females_i]['Beta']) - \
             np.array(self.results[males_i]['Beta'])
        t2 = np.sqrt(np.power(np.array(self.results[females_i]['SE']), 2) + \
                     np.power(np.array(self.results[males_i]['SE']), 2))
        zdiff = t1 / t2
        pval  = 2*stats.norm.cdf(-abs(zdiff))
        mindex   = self.results[0].index
        final_df = pd.DataFrame(index=mindex)
        final_df['SE']   = t2
        final_df['Beta'] = t1
        final_df['pvalue'] = pval
        final_df['Variable_pvalue'] = pval

        final_df.loc[~final_df['pvalue'].isna(), 'pvalue_bonferroni'] = multipletests(
                                                            final_df.loc[~final_df['pvalue'].isna(), 'pvalue'],
                                                            method='bonferroni')[1]

        self.results.append(final_df)
        self.name.append('sex differences')

    def apply_decision_tree(self):
        '''
        apply the decision tree
        '''
        #### DECISION TREE
        # 1. Z_diff is significant, also both meta female and meta male are significant, and betas are opposite
        zdiff_sign = self.results[7]['pvalue_bonferroni'] < 0.05
        both_sign  = ( self.results[4]['pvalue'] < 0.05 ) & ( self.results[5]['pvalue'] < 0.05 )
        opposite_direction = self.results[4]['Beta'] * self.results[5]['Beta'] < 0
        keep_qual  = zdiff_sign & both_sign & opposite_direction

        # 2. Overall nominal significance, zdiff significance boferroni, both significant and same direction
        overall_nominal  = self.results[6]['pvalue'] < 0.05
        zdiff_bonferroni = self.results[7]['pvalue'] < ( 0.05/ sum(overall_nominal))
        same_direction   = self.results[4]['Beta'] * self.results[5]['Beta'] > 0
        keep_quant = overall_nominal & zdiff_bonferroni & same_direction & both_sign

        # 3. Overall nominal significance, zdiff significance boferroni, only one significant
        one_sig  = ( ( self.results[4]['pvalue'] < 0.05 ) & ( self.results[5]['pvalue'] > 0.05 ) ) | \
                   ( ( self.results[4]['pvalue'] > 0.05 ) & ( self.results[5]['pvalue'] < 0.05 ) )
        keep_pure = overall_nominal & zdiff_bonferroni & one_sig

        for i in range(len(self.results)):
            self.results[i]['difference_type'] = 'None'
            self.results[i].loc[keep_qual,'difference_type']  = 'Qualitative'
            self.results[i].loc[keep_quant,'difference_type'] = 'Quantitative'
            self.results[i].loc[keep_pure,'difference_type']  = 'Pure'

    def add_variable_names(self, var_description, var_category):
        '''
        add human readable variable names, given in var_description and var_category
        '''
        for i in range(len(self.results)):
            index_variable  = self.results[i].index.get_level_values(level='Variable')
            index_phenotype = self.results[i].index.get_level_values(level='Phenotype')
            variable_name   = _add_var_from_dict(index_variable, var_description)
            self.results[i]['Variable_Name']  = variable_name
            phenotype_name  = _add_var_from_dict(index_phenotype, var_description)
            self.results[i]['Phenotype_Name'] = phenotype_name
            variable_category = _add_var_from_dict(index_variable, var_category)
            self.results[i]['Variable_Category'] = variable_category
            self.results[i]['Variable']  = index_variable
            self.results[i]['Phenotype'] = index_phenotype
            self.results[i] = self.results[i].set_index(['Variable','Variable_Name','Phenotype','Phenotype_Name'])

    def save_results(self, respath):
        '''
        Save PheEWAS results in respath
        '''
        os.chdir(respath)
        for i in range(len(self.results)):
            name = self.name[i] + '.csv'
            self.results[i].to_csv(name)


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

def get_phenotypes(print_descriptions=False, var_description=None, cleaned=False):
    '''
    Get the list of phenotypes and print their description if desired
    '''
    if cleaned is True:
        phenotype = ["LBXSCRINV","URXUCR","LBXSCR","LBXSATSI","LBXSAL","URXUMASI","URXUMA","LBXSAPSI","LBXSASSI","LBXSC3SI",
             "LBXSBU","LBXBAP","LBXCPSI","LBXCRP","LBXSCLSI","LBXSCH","LBDHDL","LBDLDL","LBXSGTSI","LBXSGB",
             "LBXGLU","LBXGH","LBXHCY","LBXSIR","LBXSLDSI","LBXMMA","LBXSOSSI","LBXSPH","LBXSKSI",	
             "LBXSNASI","LBXSTB","LBXSCA","LBXSTP","LBXSTR","LBXSUA","LBDBANO","LBXBAPCT",
             "LBDEONO","LBXEOPCT","LBXHCT","LBXHGB","LBDLYMNO","LBXMCHSI","LBXLYPCT","LBXMCVSI","LBXMPSI","LBDMONO",
             "LBXMOPCT","LBXPLTSI","LBXRBCSI","LBXRDW","LBDNENO","LBXNEPCT"] # I removed the ones that were deleted in the QC process
    elif cleaned is False:
        phenotype = ["LBXSCRINV","URXUCR","LBXSCR","LBXSATSI","LBXSAL","URXUMASI",
                 "URXUMA","LBXSAPSI","LBXSASSI","LBXSC3SI","LBXSBU","LBXBAP",
                 "LBXCPSI","LBXCRP","LBXSCLSI","LBXSCH","LBDHDL","LBDLDL",
                 "LBXFER","LBXSGTSI","LBXSGB","LBXGLU","LBXGH","LBXHCY",
                 "LBXSIR","LBXSLDSI","LBXMMA","LBXSOSSI","LBXSPH","LBXSKSI",
                 "LBXEPP","LBXSNASI","LBXTIB","LBXSTB","LBXSCA","LBXSTP",
                 "LBDPCT","LBXSTR","LBXSUA","LBDBANO","LBXBAPCT","LBDEONO",
                 "LBXEOPCT","LBXHCT","LBXHGB","LBDLYMNO","LBXMCHSI","LBXLYPCT",
                 "LBXMCVSI","LBXMPSI","LBDMONO","LBXMOPCT","LBXPLTSI",
                 "LBXRBCSI","LBXRDW","LBDNENO","LBXNEPCT","LBXIRN"]
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

def check_weights(weights):
    '''
    Check whether the weights (from load_weights function) mapped are in the
    columns of the weight dataset and return the variables that need those that
    are not there
    '''
    cols = weights[0].columns
    vals = set(weights[1].values())

    #Which survey weights are not present
    weights_absent = []
    for i in set(vals):
        if i not in cols:
            print(i + ' is not in weight dataset')
            weights_absent.append(i)

    #Which variables need those survey weights
    remove_vars = []
    for key, value in weights[1].items():
        for i in weights_absent:
            if i == value:
                remove_vars.append(key)

    return(remove_vars)

def _add_var_from_dict(index_var, var_dict):
    res = []
    for i in range(len(index_var)):
        if index_var[i] in var_dict.keys():
            res.append(var_dict[index_var[i]])
        else:
            res.append('')
    return(res)