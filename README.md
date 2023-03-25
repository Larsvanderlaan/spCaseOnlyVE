# Semiparametric inference for relative heterogeneous vaccine efficacy between strains in observational case-only studies

Semiparametric methods for inferring subgroup-specific relative vaccine efficacy in a partially vaccinated population against multiple strains of a virus. The implemented methods are designed for observational case-only studies with informative missingness in viral strain type due to vaccination status, pre-vaccination variables, and also post-vaccination factors such as viral load. Assuming that the relative strain-specific conditional vaccine efficacy has a known log-linear parametric form, we implement semiparametric asymptotically linear estimators of the parameters based on targeted (debiased) machine learning estimators for partially linear logistic regression models. 

See our manuscript for the setup, theory, and method:
https://arxiv.org/pdf/2303.11462.pdf


 

## Semiparametric inference for conditional odds ratio using partially linear regression
The method `spOR` implements a targeted maximum likelihood estimator (TMLE) for the partially linear logistic regression model. 
The function supports the ensemble learning pipeline `https://github.com/tlverse/sl3` for nuisance function estimation


## Semiparametric inference for conditional odds ratio with outcome missingness informed by post-treatment variables
The method `spORMissing` implements a targeted maximum likelihood estimator (TMLE) for the partially linear logistic regression model with outcome missingness that is informed by post-treatment variables.
The function supports the ensemble learning pipeline `https://github.com/tlverse/sl3` for nuisance function estimation.

## Data structure and target parameter
Consider the data-structure `(W, A, Z, Delta, Y)` where `W` is a vector of pre-treatment baseline covariates, `A` is a binary treatment assignment, and `Y` is a binary outcome. `Delta` is a binary indicator variable that takes the value 1 if the outcome `Y` is missing. `Z` is a post-treatment confounding variable that informs both the outcome missingness `Delta` and the outcome `Y` itself.  

The target parameter is the conditional odds ratio:
`OR(w): = {P(Y=1|A=1,W=w)/P(Y=0|A=1,W=w)} / {P(Y=1|A=0,W=w)/P(Y=0|A=0,W=w)}`,
which we model on the log scale by a linear function in w.

 
 
 
