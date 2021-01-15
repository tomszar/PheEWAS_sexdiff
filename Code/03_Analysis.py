#### LIBRARIES
import clarite
import modules

#### SET PATHS
paths  = modules.set_project_paths()

#### DATA
split_datasets  = modules.load_clean_data(paths[2])  
var_description = modules.load_var_information(paths[1])

# Categorize variables
# All unknown are continuous
for i in range(4):
    split_datasets[i].data = clarite.modify.categorize(split_datasets[i].data)
    var_unknown            = split_datasets[i].print_unknown_vars(var_description)
    split_datasets[i].data = clarite.modify.make_continuous(split_datasets[i].data, only=var_unknown)

## Survey weights
# Much information on the issue on survey weights can be found [here](https://wwwn.cdc.gov/nchs/nhanes/tutorials/module3.aspx)
weights_discovery   = modules.WeightData(modules.load_weights(paths[1]))
weights_replication = modules.WeightData(modules.load_weights(paths[1], False))

# Remove variables in survey weights that not coincide between map and data
weights_discovery.remove_weights_not_in_data()
weights_replication.remove_weights_not_in_data()  

# ## Define phenotypes and covariates
phenotypes = modules.get_phenotypes(cleaned=True)
covariates = modules.get_covariates()

# Divide the survey designs based on Nhanes Data
weights_discovery_divided   = weights_discovery.divide_by_nhanes(split_datasets[0], split_datasets[1])
weights_replication_divided = weights_replication.divide_by_nhanes(split_datasets[2], split_datasets[3])

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

#Run analysis
name_of_results = ['discovery females', 'discovery males', 'replication females', 'replication males']
total_results   = []
for i in range(4):
    res = split_datasets[i].run_phe_ewas(phenotypes, covariates, survey_designs[i])
    total_results.append(res)

final_results = modules.PheEWAS_Results(total_results, name_of_results)

final_results.meta_analyze('meta females', indices=(0,2))
final_results.meta_analyze('meta males', indices=(1,3))
final_results.meta_analyze('meta total', indices=(4,5))

final_results.estimate_sex_differences()
final_results.apply_decision_tree()

# Save files
final_results.save_results(paths[2])
