#### LIBRARIES
import analyze
import nhanes as nh

#### DATA
nhanes = analyze.NhanesClean()
nhanes.categorize_variables()
nhanes.create_survey_design()

## Run analysis
results = nhanes.run_phe_ewas()

## Run meta analysis
types = ['female',
         'male',
         'total']
for t in types:
    results.meta_analyze(type=t)

results.estimate_sex_differences()
results.apply_decision_tree()

# Add variable names for human readability
results.add_variable_names(nhanes.var_description, 
                           nhanes.var_category)

# Save files
results.save_results()
