# Statistical Learning in Preclinical Drug Proarrhythmic Assessment

This repository contains the R code to replicate the result in the paper 'Statistical Learning in Preclinical Drug Proarrhythmic Assessment'. The paper is to appear on the Journal of Biopharmaceutical Statistics. Please check our [preprint](https://arxiv.org/abs/2108.00543) for details.

## Code Structure

### 1. stem_logistic.R, stem_randomforest.R, wedge_logistic.R, wedge_randomforest.R

The model performance of ordinal logistic regression and ordinal random forest on the stem cell and wedge dataset. All metrics are point estimation.

### 2. stem_logistic_bootstrap.R, stem_radomforest_bootstrap.R, wedge_logistic_bootstrap.R, wedge_randomforest_bootstrap.R

The model performance measured by stratified bootstrap. All metrics are asymptotic distributions.


### 3. stem_importance.R, wedge_importance.R

The permutation predictor importance on the stem cell and wedge dataset.

### 4. postive_control_stem_logistic.R, postive_control_stem_randomforest.R, positive_control_wedge_logistic.R, positive_control_wedge_randomforest_paper.R

The control analysis. Model performance was measured condition on the correct prediction of control drugs.

### 5. model comparison.R

The comparison of three-category prediction accuracy among ordinal logistic regression, ordinal random forest, multinomial logistic regression, multinomial random forest, and classical ordinal logistic regression.

### 6. figure_paper.R

All the codes that draw figures in the paper.

### 7. Result

All the numerical results output by the code. These results generate the figures and tables in the paper.

# Contact
If you have any suggestions or comments on this study, please contact Nan Miles Xi (<mxi1@luc.edu>). 

# Citation
If you use the R code in your work, please cite 

Xi, M.N., Hsu, Y., Dang, Q., and Huang, D. (2022). Statistical Learning in Preclinical Drug Proarrhythmic Assessment. arXiv:2108.00543.
