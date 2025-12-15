# imputation_pipeline

R code implementing a class and simulation study for the imputation pipeline in my cardiac AKI project.

The code in this repository is responsible for the results in the STS 2026 abstract. Although not in this repository, there are two csv files used for this analysis overall. The first dataset titled ``primary_data_unscaled_feature_learn.csv`` is used to learn the relevant features for AKI. The second dataset titled ``primary_dataset_unscaled.csv`` is used for computing causal effects and mediation effects.

Something I still need to do is to refactor the code for learning features to be in a similar format as this where I can run a .R file and be done with it.

At a high level, there are two datasets, one larger dataset without urinary biomarker data and one smaller dataset with urinary biomarker data. We believe that the selection of patients into the two datasets is completely at random, so we can learn the distribution of the urinary biomarker data given the treatments and covariates from the smaller dataset and use that learned distribution to impute the biomarker data in the larger dataset, allowing us to utilize the treatment, covariate, and outcome data from the larger dataset as well. We want to build a class in R for performing this analysis in an easy way on the real data.

There are four urinary biomarkers. If we assume that these four biomarkers are direct measurements of kidney function (i.e. they are the M's directly), then we can use all four of these variables in the mediation analysis. 

The steps of the analysis are outlined below:

1. Learn the densities f(M_1 | M_2, M_3, M_4, A_1, A_2, X), f(M_2 | M_3, M_4, A_1, A_2, X), f(M_3 | M_4, A_1, A_2, X), f(M_4 | A_1, A_2, X) from the small dataset.
2. Using the learned densities in the previous step, sample/impute biomarker values M_1-M_4 in the larger dataset. All of the following steps now use both datasets and the maximum amount of data.
3. Learn the densities f(A_1 | A_2, X), f(A_2 | X), f(A_1 | A_2), and f(A_2).
4. Compute MSM weights f(A_1 | A_2) f(A_2) / f(A_1 | A_2, X) f(A_2 | X). The values in the numerator are the observed marginals.
5. Compute the pseudo-outcomes Y' = Y f(M_1, M_2, M_3, M_4 | A_1=a_1, A_2=a_2, X) / f(M_1, M_2, M_3, M_4 | A_1=a_1', A_2=a_2', X).
6. Fit a model for the pseudo-outcome Y' given A_1 and A_2 using the MSM weights and evaluate it at (A_1, A_2) = (a_1', a_2'). This gives an estimate for the mixed mediation term.
7. Fit a model for the outcome Y given A_1 and A_2 using the MSM weights and evaluate it at (A_1, A_2) = (a_1, a_2). This gives an estimate for the mean of the counterfactual Y(a_1, a_2).
8. Take the difference between steps 6 and 7 for the direct effect of oxygen delivery on AKI.

2025-11-19
Additional analyses
1. Do the same analysis as before but only use one biomarker at a time and see what happens.
2. Take the mean of the delta change in biomarkers at all 4 time points and repeat the same analysis as before.

These two additional analyses should not be too difficult. The first one involves changing the data dictionary in the data_application.R file. We will create four new files that are copies of the original ``data_application.R`` where each file only uses one mediator.
The second involves computing new columns for the biomarker data and then pushing through the new dataset through the same code. After computing the new column, we will create a copy of ``data_application.R`` that uses data for these new columns as the mediator.

2025-12-03
Additional analyses
1. Compute the causal effect of causally significant treatments on each of the four biomarkers.

This analysis can just use the smaller dataset with biomarker data. This is because using an imputation method will not increase the efficiency of the method. One additional task is that we will need to implement the computing the counterfactual mean for a continuous outcome. Actually, it is already implemented, so hopefully it won't be too complicated.

2025-12-15
Note to self:
- Note that the ``data_application.R`` file in this repo relies on data from a file titled ``primary_data_final.csv``. This file is obtained from ``cardiac_surgery_aki/feature_selection/feature_selection.Rmd`` where I performed a sample split. The first sample split was used to do feature selection on which of the 17 selected features are potentially causally relevant. The second sample split was used to estimate causal effects.
- The ``feature_selection.Rmd`` file then further depends on ``perfusion_data.csv``, ``combined_data.csv``, and ``DO2_curve_features.csv``. The ``perfusion_data.csv`` ``combined_data.csv`` files come from ``datamerge/data_merge.Rmd`` where I merge the perfusion and flatfile together. The ``DO2_curve_features.csv`` file comes from ``data_exploration/functional_data_analysis.Rmd`` where I compute the 17 features that we think are plausible causes of AKI. Detailed updates are from the February 2, 2025 commit of ``cardiac_surgery_aki`` repo.
