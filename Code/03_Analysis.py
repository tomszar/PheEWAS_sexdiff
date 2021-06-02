#### LIBRARIES
import clarite
import modules

#### SET PATHS
paths  = modules.set_project_paths()

#### DATA
split_datasets  = modules.load_clean_data(paths[2])  
var_description = modules.load_var_information(paths[1])
var_category    = modules.load_var_information(paths[1], description=False)

# Categorize variables
# All unknown are continuous
for i in range(4):
    split_datasets[i].data = clarite.modify.\
                             categorize(split_datasets[i].data)
    var_unknown            = split_datasets[i].\
                             print_unknown_vars(var_description)
    split_datasets[i].data = clarite.modify.\
                             make_continuous(split_datasets[i].data,
                             only=var_unknown)

## Survey weights
# Much information on the issue on survey weights can be found 
# [here](https://wwwn.cdc.gov/nchs/nhanes/tutorials/module3.aspx)
weights_discovery   = modules.WeightData(modules.load_weights(paths[1]))
weights_replication = modules.WeightData(modules.load_weights(paths[1], False))

# Remove variables in survey weights that not coincide between map and data
weights_discovery.remove_weights_not_in_data()
weights_replication.remove_weights_not_in_data()  

# ## Define phenotypes and covariates
phenotypes = modules.get_phenotypes(cleaned=True)
covariates = modules.get_covariates()

# Divide the survey designs based on Nhanes Data
weights_discovery_divided = weights_discovery.\
                            divide_by_nhanes(split_datasets[0],
                            split_datasets[1])
weights_replication_divided = weights_replication.\
                              divide_by_nhanes(split_datasets[2],
                              split_datasets[3])

survey_designs = []
for i in range(4):
    if i < 2:
        ind = i
        survey = weights_discovery_divided[ind].create_survey_design()
        survey_designs.append(survey)
    else:
        ind = i - 2
        survey = weights_replication_divided[ind].create_survey_design()
        survey_designs.append(survey)

# Run analysis
total_results = []
for i in range(4):
    res = split_datasets[i].run_phe_ewas(phenotypes,
                                         covariates, 
                                         survey_designs[i])
    total_results.append(res)

order_files = ['Discovery_Females',
               'Discovery_Males',
               'Replication_Females',
               'Replication_Males']

for i in range(len(total_results)):
    n_converged     = sum(total_results[i]['Converged'])
    n_not_converged = sum(total_results[i]['Converged']==False)
    print('For the ' + 
          order_files[i] + 
          ' dataset, there are ' +
          str(n_converged) + 
          ' converged tests, and ' +
          str(n_not_converged) + 
          ' non-converged\n')

final_results = modules.PheEWAS_Results(total_results)

final_results.meta_analyze(type='female')
final_results.meta_analyze(type='male')
final_results.meta_analyze(type='total')

final_results.estimate_sex_differences()
final_results.apply_decision_tree()

# Add variable names for human readability
final_results.add_variable_names(var_description, var_category)

# Save files
final_results.save_results(paths[2])
