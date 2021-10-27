#### IMPORT MODULES
import clean

#### READ DATA
nhanes = clean.NhanesRaw()

#### QC PROCESS
nhanes.categorize_variables()
nhanes.clarite_filters()
nhanes.transform_continuous_log2()
nhanes.zscore_normalization()
nhanes.remove_continuous_outliers()

#Run categorize and filters again because NAs 
#were entered in removal of continuous outliers
nhanes.categorize_variables()
nhanes.clarite_filters()

#Search for perfect correlation between variables
nhanes.remove_perfect_correlations()

#Plot and save
nhanes.plot_variables()
nhanes.save_data()
