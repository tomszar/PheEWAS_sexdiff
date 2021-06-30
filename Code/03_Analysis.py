#### LIBRARIES
import clarite
import nhanes as nh

#### DATA
paths   = nh.set_project_paths()
nhanes  = nh.NhanesData(paths[1], paths[2], raw=False)
weights = nh.WeightData(paths[1])
weights.harmonize()
nhanes.print_how_many_types()

## Categorize variables, make unkown continuous, and force
# phenotypes as continuous
nhanes.categorize_variables()
nhanes.make_unknown_continuous()
nhanes.make_phenotypes_continuous()

## Harmonize weights and create survey designs
weights.harmonize_with_nhanes(nhanes)
weights.create_survey_design()

## Run analysis
results = nhanes.run_phe_ewas(weights)

## Run meta analysis
results.meta_analyze(type='female')
results.meta_analyze(type='male')
results.meta_analyze(type='total')

results.estimate_sex_differences()
results.apply_decision_tree()

# Add variable names for human readability
results.add_variable_names(nhanes.var_description, 
                           nhanes.var_category)

# Save files
results.save_results(paths[2])
