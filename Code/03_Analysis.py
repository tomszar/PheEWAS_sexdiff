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
nhanes_females = pd.read_csv('CleanData_Females.csv').set_index('ID')
nhanes_males = pd.read_csv('CleanData_Males.csv').set_index('ID')
os.chdir(os.path.join(datapath, 'nh_99-06'))
var_description = pd.read_csv('VarDescription.csv')                     .drop_duplicates()                     .set_index('var')
# Convert variable descriptions to a dictionary for convenience
var_descr_dict = var_description['var_desc'].to_dict()


# ## Categorize variables

# In[4]:


nhanes_females = cl.modify.categorize(nhanes_females)
nhanes_males   = cl.modify.categorize(nhanes_males)


# In[5]:


var_types_f   = cl.describe.get_types(nhanes_females)
var_unknown_f = var_types_f[var_types_f == 'unknown'].index

var_types_m   = cl.describe.get_types(nhanes_males)
var_unknown_m = var_types_m[var_types_m == 'unknown'].index

print('Unknown variables in females:')
for v in list(var_unknown_f):
    if v in var_descr_dict:
        print(f"\t{v} = {var_descr_dict[v]}")
        
print('\nUnknown variables in males:')
for v in list(var_unknown_m):
    if v in var_descr_dict:
        print(f"\t{v} = {var_descr_dict[v]}")


# In[6]:


nhanes_females = cl.modify.make_continuous(nhanes_females, only=var_unknown_f)
nhanes_males   = cl.modify.make_continuous(nhanes_males, only=var_unknown_m)


# ## Survey weights
# 
# Much information on the issue on survey weights can be found [here](https://wwwn.cdc.gov/nchs/nhanes/tutorials/module3.aspx)

# In[7]:


os.chdir(datapath)
survey_design = pd.read_csv("weights_8yr.csv", sep=",")                            .rename(columns={'SEQN':'ID'})                            .set_index("ID")
survey_design.head()


# In[8]:


os.chdir(datapath)
# These files map variables to their correct weights, and were compiled by reading throught the NHANES codebook
var_weights = pd.read_csv("VarWeights8yr.csv")
var_weights.head()


# In[9]:


# Convert the data to two dictionaries for convenience
weights = var_weights.set_index('variable_name')['weight'].to_dict()


# ## Define phenotypes and covariates
# 

# In[10]:


phenotype = ["LBXSCRINV","URXUCR","LBXSCR","LBXSATSI","LBXSAL","URXUMASI","URXUMA","LBXSAPSI","LBXSASSI","LBXSC3SI",
             "LBXSBU","LBXBAP","LBXCPSI","LBXCRP","LBXSCLSI","LBXSCH","LBDHDL","LBDLDL","LBXFER","LBXSGTSI","LBXSGB",
             "LBXGLU","LBXGH","LBXHCY","LBXSIR","LBXSLDSI","LBXMMA","LBXSOSSI","LBXSPH","LBXSKSI","LBXEPP",	
             "LBXSNASI","LBXTIB","LBXSTB","LBXSCA","LBXSTP","LBDPCT","LBXSTR","LBXSUA","LBDBANO","LBXBAPCT",
             "LBDEONO","LBXEOPCT","LBXHCT","LBXHGB","LBDLYMNO","LBXMCHSI","LBXLYPCT","LBXMCVSI","LBXMPSI","LBDMONO",
             "LBXMOPCT","LBXPLTSI","LBXRBCSI","LBXRDW","LBDNENO","LBXNEPCT","LBXIRN"]
for v in phenotype:
    if v in var_descr_dict:
        print(f"\t{v} = {var_descr_dict[v]}")
        
covariates = ["black", "mexican", "other_hispanic", "other_eth", "SES_LEVEL", "RIDAGEYR", "SDDSRVYR","BMXBMI"]


# In[11]:


#### REMOVE DROPED PARTICIPANTS FROM NHANES TO SURVEY
survey_design_females = survey_design.loc[nhanes_females.index]
survey_design_males    = survey_design.loc[nhanes_males.index]

sd_desingspec_females = cl.survey.SurveyDesignSpec(survey_df=survey_design_females,
                                                   strata="SDMVSTRA",
                                                   cluster="SDMVPSU",
                                                   nest=True,
                                                   weights=weights,
                                                   single_cluster='adjust',
                                                   drop_unweighted= True,
                                                   fpc=None)

sd_desingspec_males = cl.survey.SurveyDesignSpec(survey_df=survey_design_males,
                                                  strata="SDMVSTRA",
                                                  cluster="SDMVPSU",
                                                  nest=True,
                                                  weights=weights,
                                                  single_cluster='adjust',
                                                  drop_unweighted= True,
                                                  fpc=None)


# In[12]:


results_females = []
results_males   = []
for current_pheno in phenotype:
    print(current_pheno)
    drop_phenos = [p for p in phenotype if p != current_pheno]
    ewas_result_dis = cl.analyze.ewas(phenotype  = current_pheno,
                                      covariates = covariates,
                                      min_n = 200,
                                      data = nhanes_females.drop(columns=drop_phenos),
                                      regression_kind = 'weighted_glm',
                                      survey_design_spec=sd_desingspec_females)
    results_females.append(ewas_result_dis)
    ewas_result_dis = cl.analyze.ewas(phenotype  = current_pheno,
                                      covariates = covariates,
                                      min_n = 200,
                                      data = nhanes_males.drop(columns=drop_phenos),
                                      regression_kind = 'weighted_glm',
                                      survey_design_spec=sd_desingspec_males)
    results_males.append(ewas_result_dis)
    
results_females = pd.concat(results_females)
results_males   = pd.concat(results_males)
cl.analyze.add_corrected_pvalues(results_females)
cl.analyze.add_corrected_pvalues(results_males)


# ## Sex difference analysis

# In[74]:


# Concatenate columns
results = pd.merge(results_females, results_males, how='inner', left_index=True, right_index=True, suffixes=('_female', '_male'))

# Kepp only those converged in both males and females
keep_results = (results.Converged_female & results.Converged_male)
results      = results.loc[keep_results]

# Estimate differences
t1 = np.array(results['Beta_female']) - np.array(results['Beta_male'])
t2 = np.sqrt(np.power(np.array(results['SE_female']), 2) + np.power(np.array(results['SE_male']), 2))
zdiff   = t1 / t2
pvalues = 2*stats.norm.cdf(-abs(zdiff))
results['pvalue_difference'] = pvalues

# Save file
os.chdir(respath)
results.to_csv('Difference_test.csv')


# In[ ]:




