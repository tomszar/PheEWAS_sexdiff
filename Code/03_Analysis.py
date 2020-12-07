#!/usr/bin/env python
# coding: utf-8

# # Analysis
# 
# ## Loading data

# In[1]:


import os
import pandas as pd
import clarite as cl
import numpy as np
import scipy.stats as stats


# In[2]:


#### SET PATHS
os.chdir('..')
mainpath = os.getcwd()
datapath = os.path.join(mainpath, 'Data')
respath  = os.path.join(mainpath, 'Results')


# In[3]:


os.chdir(respath)
loadfiles = ['Discovery_Females', 'Discovery_Males', 'Replication_Females', 'Replication_Males']
dfs = []
for i in range(4):
    dfs.append(pd.read_csv('CleanData_' + loadfiles[i] + '.csv').set_index('ID'))

os.chdir(os.path.join(datapath, 'nh_99-06'))
var_description = pd.read_csv('VarDescription.csv')                     .drop_duplicates()                     .set_index('var')

# Convert variable descriptions to a dictionary for convenience
var_descr_dict = var_description['var_desc'].to_dict()


# ## Categorize variables

# In[4]:


for i in range(len(dfs)):
    dfs[i] = cl.modify.categorize(dfs[i])


# In[5]:


text  = ['Discovery females', 'Discovery males', 'Replication females', 'Replication males']
for i in range(len(dfs)):
    var_types   = cl.describe.get_types(dfs[i])
    var_unknown = var_types[var_types == 'unknown'].index
    print('Unknown variables in ' + text[i] + ':')
    for v in list(var_unknown):
        if v in var_descr_dict:
            print(f"\t{v} = {var_descr_dict[v]}")


# In[6]:


for i in range(len(dfs)):
    var_types   = cl.describe.get_types(dfs[i])
    var_unknown = var_types[var_types == 'unknown'].index
    dfs[i] = cl.modify.make_continuous(dfs[i], only=var_unknown)


# In[7]:


for i in range(len(dfs)):
    cl.describe.summarize(dfs[i])


# ## Survey weights
# 
# Much information on the issue on survey weights can be found [here](https://wwwn.cdc.gov/nchs/nhanes/tutorials/module3.aspx)

# In[8]:


os.chdir(datapath)
survey_design_discovery = pd.read_csv("weights_discovery.txt", sep="\t")                            .rename(columns={'SEQN':'ID'})                            .set_index("ID")
survey_design_discovery.head()


# In[9]:


os.chdir(datapath)
survey_design_replication = pd.read_csv("weights_replication.txt", sep="\t")                            .rename(columns={'SEQN':'ID'})                            .set_index("ID")
survey_design_replication.head()


# In[10]:


# Divide replication weights by 2 to get 4 year weights
survey_design_replication.iloc[:,3:] = survey_design_replication.iloc[:,3:] / 2
survey_design_replication.head()


# In[11]:


os.chdir(datapath)
# These files map variables to their correct weights, and were compiled by reading throught the NHANES codebook
var_weights = pd.read_csv("VarWeights.csv")
var_weights.head()


# In[12]:


# Convert the data to two dictionaries for convenience
weights_discovery   = var_weights.set_index('variable_name')['discovery'].to_dict()
weights_replication = var_weights.set_index('variable_name')['replication'].to_dict()

# Remove the list of variables from replication (trows an error and it's not in the dataset either)
#remove_rep_dict = ['LBXTSH', 'URX24D', 'URX25T', 'URX4FP', 'URXACE', 'URXATZ', 'URXCB3', 'URXCCC', 'URXCMH', 'URXCPM',
#                   'URXDEE', 'URXDIZ', 'URXDPY', 'URXMET', 'URXOPM']
#for i in remove_rep_dict:
#    weights_replication.pop(i, None)


# In[13]:


cols = survey_design_replication.columns
vals = set(weights_replication.values())

#Which survey weights are not present
weights_absent = []
for i in set(vals):
    if i not in cols:
        print(i + ' is not in cols')
        weights_absent.append(i)

#Which variables need those survey weights
remove_vars = []
for key, value in weights_replication.items():
    for i in weights_absent:
        if i == value: 
            remove_vars.append(key)
            
# Remove those vars from the weights dict
for i in remove_vars:
    weights_discovery.pop(i)
    weights_replication.pop(i)


# ## Define phenotypes and covariates
# 

# In[14]:


phenotype = ["LBXSCRINV","URXUCR","LBXSCR","LBXSATSI","LBXSAL","URXUMASI","URXUMA","LBXSAPSI","LBXSASSI","LBXSC3SI",
             "LBXSBU","LBXBAP","LBXCPSI","LBXCRP","LBXSCLSI","LBXSCH","LBDHDL","LBDLDL","LBXSGTSI","LBXSGB",
             "LBXGLU","LBXGH","LBXHCY","LBXSIR","LBXSLDSI","LBXMMA","LBXSOSSI","LBXSPH","LBXSKSI",	
             "LBXSNASI","LBXSTB","LBXSCA","LBXSTP","LBXSTR","LBXSUA","LBDBANO","LBXBAPCT",
             "LBDEONO","LBXEOPCT","LBXHCT","LBXHGB","LBDLYMNO","LBXMCHSI","LBXLYPCT","LBXMCVSI","LBXMPSI","LBDMONO",
             "LBXMOPCT","LBXPLTSI","LBXRBCSI","LBXRDW","LBDNENO","LBXNEPCT"] # I removed the ones that were deleted in the QC process
for v in phenotype:
    if v in var_descr_dict:
        print(f"\t{v} = {var_descr_dict[v]}")
        
covariates = ["black", "mexican", "other_hispanic", "other_eth", "SES_LEVEL", "RIDAGEYR", "SDDSRVYR","BMXBMI"]


# In[15]:


# Create 4 survey designs
surveys = [survey_design_discovery.loc[dfs[0].index], survey_design_discovery.loc[dfs[1].index],
          survey_design_replication.loc[dfs[2].index], survey_design_replication.loc[dfs[3].index]]
sds = []
for i in range(len(dfs)):
    if i == 0 or i == 1:
        w = weights_discovery
    else:
        w = weights_replication
        
    sds.append(cl.survey.SurveyDesignSpec(survey_df=surveys[i],
                                                   strata="SDMVSTRA",
                                                   cluster="SDMVPSU",
                                                   nest=True,
                                                   weights=w, # This are the same in discovery and replication
                                                   single_cluster='adjust'))


# In[16]:


#Omit pandas warning
import warnings
warnings.filterwarnings("ignore")

#Run analysis
results_d_f = []
results_d_m = []
results_r_f = []
results_r_m = []
for i in range(len(dfs)):
    for current_pheno in phenotype:
        print(current_pheno)
        drop_phenos = [p for p in phenotype if p != current_pheno]
        ewas_result_dis = cl.analyze.ewas(phenotype = current_pheno,
                                          covariates = covariates,
                                          min_n = 200,
                                          data = dfs[i].drop(columns=drop_phenos),
                                          regression_kind = 'weighted_glm',
                                          survey_design_spec=sds[i])
        if i == 0:
            results_d_f.append(ewas_result_dis)
        elif i == 1:
            results_d_m.append(ewas_result_dis)
        elif i == 2:
            results_r_f.append(ewas_result_dis)
        elif i == 3:
            results_r_m.append(ewas_result_dis)
                
    
results = [pd.concat(results_d_f).merge(pd.concat(results_r_f), how = 'inner', 
                             left_index=True, right_index=True, suffixes =['_d', '_r'] ), 
           pd.concat(results_d_m).merge(pd.concat(results_r_m), how = 'inner', 
                             left_index=True, right_index=True, suffixes =['_d', '_r'] )]


# ## Meta-analyze

# In[17]:


def meta_analysis(betas, ses):
    '''
    Perform a meta-analysis from a given set of beta and se coefficients
    '''
    import numpy as np
    import scipy.stats as stats
    
    w1 = 1 / np.power(ses[0],2)
    w2 = 1 / np.power(ses[1],2)
    meta_se   = np.sqrt( 1 / (w1 + w2))
    meta_beta = (betas[0] * w1 + betas[1] * w2) / (w1 + w2)
    zeta = meta_beta / meta_se
    pval = 2 * stats.norm.cdf(-abs(zeta))
    
    return[meta_se, meta_beta, pval]


# In[18]:


# Males vs females meta-analysis
drop_cols = ['Variable_type_r', 'Weight_r', 'Converged_r', 'LRT_pvalue_r', 'Diff_AIC_r'] #Columns to drop

for i in range(len(results)): # Calculate the meta-analysis parameters
    vals = meta_analysis([results[i]['Beta_d'], results[i]['Beta_r']] , [results[i]['SE_d'], results[i]['SE_r']])
    results[i]['SE']     = vals[0]
    results[i]['Beta']   = vals[1]
    results[i]['pvalue'] = vals[2]
    results[i]['Variable_pvalue'] = vals[2]
    results[i]['N'] = results[i]['N_d'] + results[i]['N_r']
    results[i] = results[i].rename(columns={'Variable_type_d':'Variable_type',
                                            'Weight_d': 'Weight',
                                            'Converged_d':'Converged',
                                            'LRT_pvalue_d':'LRT_pvalue',
                                            'Diff_AIC_d':'Diff_AIC'})
    # Remove some unnecesary columns
    results[i] = results[i].drop(drop_cols, axis=1)


# In[19]:


# Final meta analysis (males with females)
vals = meta_analysis([results[0]['Beta'], results[1]['Beta'] ], [results[0]['SE'], results[1]['SE']] )

final_meta = pd.DataFrame(results[0].iloc[:,0:3], index=results[0].index)
final_meta['N'] = results[0]['N'] + results[1]['N']
final_meta['LRT_pvalue'] = results[0]['LRT_pvalue']
final_meta['Diff_AIC']   = results[0]['Diff_AIC']
final_meta['SE']     = vals[0]
final_meta['Beta']   = vals[1]
final_meta['pvalue'] = vals[2]
final_meta['Variable_pvalue'] = vals[2]
cl.analyze.add_corrected_pvalues(final_meta)


# ## Sex difference analysis
# 
# There are two complementary analysis done here:
# - First, a wide difference test (better suited to detect opposite effects), that will take all possible tests, calculate bonferroni corrected pvalues, and then only keep those that have opposite effects between sexes
# - Second, a filtering first test (better suited to detect effects that are only or more pronounced in one sex), that will take the overall meta-analyzed association, filtering by (what threshold?), and then test for differences between sexes and retain only those that are not opposite effects

# In[20]:


# Merge datasets and rename and remove columns
drop_cols   = ['Variable_type_m', 'Weight_m', 'Converged_m', 'LRT_pvalue_m', 'Diff_AIC_m',
               'Variable_type_f', 'Weight_f', 'Converged_f', 'LRT_pvalue_f', 'Diff_AIC_f'] #Columns to drop

merge_sex = pd.merge(results[0], results[1], how = 'inner', left_index=True, right_index=True, suffixes =['_f', '_m'])

# Kepp only those converged in both males and females
keep_results = (merge_sex['Converged_f'] & merge_sex['Converged_m'])
merge_sex    = merge_sex.loc[keep_results]

# Remove and rename columns
merge_sex = merge_sex.drop(drop_cols, axis=1)

# Merge with overall meta analysis
final_result = pd.merge(final_meta, merge_sex, how = 'inner', left_index=True, right_index=True, suffixes=[None, 'y'])


# In[36]:


# Estimate differences
t1 = np.array(final_result['Beta_f']) - np.array(final_result['Beta_m'])
t2 = np.sqrt(np.power(np.array(final_result['SE_f']), 2) + np.power(np.array(final_result['SE_m']), 2))
zdiff   = t1 / t2
pval = 2*stats.norm.cdf(-abs(zdiff))
final_result['pvalue_diff'] = pval
final_result['SE_diff'] = t2
final_result['Beta_diff'] = t1
final_result['Variable_pvalue_diff'] = pval

# Get corrected pvalues
from statsmodels.stats.multitest import multipletests
final_result.loc[~final_result['pvalue_diff'].isna(), 'pvalue_diff_bonferroni'] = multipletests(
                                                                                        final_result.loc[~final_result['pvalue_diff'].isna(), 'pvalue_diff'], 
                                                                                        method='bonferroni'
                                                                                        )[1]

# Save file
os.chdir(respath)
final_result.to_csv('Difference_test.csv')

