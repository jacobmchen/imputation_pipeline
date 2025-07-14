# imputation_pipeline

R code implementing a class and simulation study for the imputation pipeline in my cardiac AKI project.

At a high level, there are two datasets, one larger dataset without urinary biomarker data and one smaller dataset with urinary biomarker data. We believe that the selection of patients into the two datasets is completely at random, so we can learn the distribution of the urinary biomarker data given the treatments and covariates from the smaller dataset and use that learned distribution to impute the biomarker data in the larger dataset, allowing us to utilize the treatment, covariate, and outcome data from the larger dataset as well. We want to build a class in R for performing this analysis in an easy way on the real data.

There are four urinary biomarkers. If we assume that these four biomarkers are direct measurements of kidney function (i.e. they are the M's directly), then we can use all four of these variables in the mediation analysis. 

The steps of the analysis are outlined below:

1. Learn the densities f(M_1 | M_2, M_3, M_4, A_1, A_2, X), f(M_2 | M_3, M_4, A_1, A_2, X), f(M_3 | M_4, A_1, A_2, X), f(M_4 | A_1, A_2, X) from the small dataset.
2. Using the learned densities in the previous step, sample/impute biomarker values M_1-M_4 in the larger dataset. All of the following steps now use both datasets and the maximum amount of data.
3. Learn the densities f(A_1 | A_2, X), f(A_2 | X), f(A_1 | A_2), and f(A_2).
4. Compute MSM weights f(A_1 | A_2) f(A_2) / f(A_1 | A_2, X) f(A_2 | X).
5. Compute the pseudo-outcomes Y' = Y f(M_1, M_2, M_3, M_4 | A_1=a_1, A_2=a_2, X) / f(M_1, M_2, M_3, M_4 | A_1=a_1', A_2=a_2', X).
6. Fit a model for the pseudo-outcome Y' given A_1 and A_2 using the MSM weights and evaluate it at (A_1, A_2) = (a_1', a_2'). This gives an estimate for the mixed mediation term.
7. Fit a model for the outcome Y given A_1 and A_2 using the MSM weights and evaluate it at (A_1, A_2) = (a_1, a_2). This gives an estimate for the mean of the counterfactual Y(a_1, a_2).
8. Take the difference between steps 6 and 7 for the direct effect of oxygen delivery on AKI.