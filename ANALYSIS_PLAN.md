# Sex differences in PheEWAS analysis plan

The purpose of this analysis is to discover phenome-environment-wide associations that have differential effects between sexes. Most of the pipeline, phenotype and variable selection will follow [Nikki’s plan](https://docs.google.com/document/d/1_2FWZHSnPEc1CqDxdVDRocUD1A3srALysvxBdCpunqY/edit?usp=sharing)

1. Initial QC process
    - Drop any variables that are indeterminate according to the NHANES data dictionary
    - Remove variables that don’t have 4 year weights
    - Drop any non-environmental exposures (examples physical fitness)
    - Determine covariates and phenotypes
    - From the remainder of the variables split them based on their type
    - Manually inspect the definition of ambiguous variables and determine their type
    - Merge binary into categorical variables
2. Make an adjustment to the lipid variable (LBDLDL) based on statin medication
3. Split by adults (greater than 18 yo) and subadults, and use only adults thereafter
4. Split dataset into discovery (series 1 and 2) and replication (series 3 and 4) datasets
    - Make sure the phenotypes from the dataset categorized under biochemistry and blood lab test/measures are included in at least one series from Discovery and Replication
    - Nikki’s list of 58 phenotypes already has taken care of the previous step
5. Split datasets into males and females (therefore, ending up with 4 different data sets)
6. Drop any variable that has a missing value in the covariate list
7. Separate datasets into phenotypes and exposures and continue with QC
8. Remove phenotypes with more than 90% of samples with 0 value 
    - Nikki’s list of phenotypes already has those removed
9. Keep only those phenotypes and variables that are present in all 4 data sets
10. Normalize variables by mean subtraction and sd division

Given the complexity of the EWAS models, it is easier and more convenient to run a stratified EWAS and test for sex differences after. Therefore, we will run four separate EWAS models for each dataset independently and estimating the corresponding parameters ($\beta$, $se$). Winkler et al (2017) recommend following two approaches in parallel in genome-wide association studies if there is no prior hypothesis on sex differences: to run a genome-wide difference test between sexes to search for opposite effects, and another approach that first filters for an overall association and then test for the difference between sexes to search for those with differences in the size of the effect, or for those that there is no effect in one sex. The following categories will be used to refer to those differences in effects:

- Qualitative: exposures that have opposite effects between sexes
- Quantitative: exposures have the same direction of effect, but the effect is larger in one sex
- Pure: exposures have an effect in only one sex and not in the other

## Phenome-environment-wide sex difference test

Considering two sexes, $i=1,2$, let $Y_p$ be a vector of phenotypes, where $p=1...,P$ considering $P$ phenotypes, and $X_q$ is a vector of environmental exposures, $q=1...,Q$, considering $Q$ environmental exposures. We write the linear regression as:

$$
Y_{ip} = X_{iq}\beta + \epsilon
$$ 

Our interest will be focused on $\beta_{ipq}$ which is the beta coefficient of the effect of the environmental exposure $q$ on phenotype $p$, in sex $i$, with its corresponding standard error $se_{ipq}$. For the phenome-environment-wide sex difference test we will estimate the *difference test* as:

$$
Z_{diff} = \frac{\beta_{1pq} - \beta_{2pq}}{\sqrt{se^2_{1pq} + se^2_{2pq}}}
$$

We will use the Bonferroni correction for multiple testing.

## Filtering by overall association

The second pipeline we will use to estimate sex differences incorporates a filtering based on overall association before the difference test. On a stratified approach, the *overall test* is given by:

$$
Z_{overall} = \frac{\frac{\beta_{1pq} }{se^2_{1pq} } + \frac{\beta_{2pq} }{se^2_{2pq} } }{\sqrt{ \frac{1}{se^2_{1pq} } + \frac{1}{se^2_{2pq} } } }
$$

Both the filtering by overall test and the subsequent difference test will use the Bonferroni correction for multiple testing.

## Step-by-step analysis

Finally, the pipeline will follow these steps:

1. To detect qualitative effect differences we will run the difference test and select those with a Bonferroni corrected $\alpha$. Then, we will keep only those with opposite effects (different directions of $\beta$ coefficients, and nominally significant in both sexes)
2. To detect quantitative and pure effect differences we will first filter by an overall effect $p-value < 0.05$ and then select those with a difference test with a Bonferroni corrected $\alpha$. We will classify those will $\beta$ coefficients in the same direction and with a nominal p-value in both sexes as a quantitative difference, while those with a nominal p-value in only one sex will be classified as a pure difference

# References

Winkler, Thomas W., Anne E. Justice, L. Adrienne Cupples, Florian Kronenberg, Zoltán Kutalik, Iris M. Heid, and GIANT consortium. 2017. “Approaches to Detect Genetic Effects That Differ between Two Strata in Genome-Wide Meta-Analyses: Recommendations Based on a Systematic Evaluation.” PloS One 12 (7): e0181038. https://doi.org/10/gbqbnm.
