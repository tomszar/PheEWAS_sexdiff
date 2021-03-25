#### LIBRARIES
import clarite
import modules
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#### SET PATHS
paths  = modules.set_project_paths()

#### READ DATA
final_results = modules.load_results(paths[2])
clean_data    = modules.load_clean_data(paths[2])
var_description = modules.load_var_information(paths[1])
var_category    = modules.load_var_information(paths[1], description=False)
phenotypes      = modules.get_phenotypes(cleaned=True)

#### ADD CLARITE REQUIREMENTS
for i in (4,5,7):
    final_results.results[i]['Converged']     = True
    final_results.results[i]['Converged'][final_results.results[i]['pvalue'].isnull()] = False
    if i == 4:
        pval_t = 1e-299
    elif i == 5:
        pval_t = 1e-303
    elif i == 7:
        final_results.results[i]['N'] = np.array(final_results.results[4]['N']) + np.array(final_results.results[5]['N'])
    final_results.results[i].loc[final_results.results[i]['pvalue'] == 0, 'pvalue'] = pval_t
    final_results.results[i]['Variable_type'] = final_results.results[0]['Variable_type']
    final_results.results[i]['Diff_AIC']   = final_results.results[0]['Diff_AIC']
    final_results.results[i]['LRT_pvalue'] = final_results.results[0]['LRT_pvalue']
    final_results.results[i]['pvalue_bonferroni'] = final_results.results[i]['pvalue']
    final_results.results[i]['pvalue_fdr']        = final_results.results[i]['pvalue']
    final_results.results[i] = final_results.results[i].set_index(['Variable','Phenotype'])

#### MANHATTAN PLOTS
os.chdir(paths[2])
fig = plt.figure(figsize=(12,8), dpi=300)
clarite.plot.manhattan({'Females':final_results.results[4], 'Males':final_results.results[5]}, 
                                    num_labeled=0,
                                    categories=var_category, 
                                    title="Weighted EWAS Sex Difference Results",
                                    figure=fig,
                                    filename='Plots/FigS1_ManhattanMalesFemales.png')

fig = plt.figure(figsize=(12,6), dpi=300)
clarite.plot.manhattan({'Total':final_results.results[7]}, 
                                    num_labeled=0,
                                    categories=var_category, 
                                    title="Weighted EWAS Sex Difference Results",
                                    figure=fig,
                                    filename='Plots/Fig1_ManhattanSexDifferences.png',
                                    bonferroni=0.1)

