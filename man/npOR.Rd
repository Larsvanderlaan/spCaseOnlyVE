% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npOR.R
\name{npOR}
\alias{npOR}
\title{Nonparametric Targeted inference for the conditional odds ratio with outcome missingness.}
\usage{
npOR(
  working_formula = logOR ~ 1,
  W,
  A,
  Y,
  Delta = NULL,
  weights = NULL,
  W_new = W,
  glm_formula_Y0W = NULL,
  sl3_learner_Y0W = NULL,
  glm_formula_OR = NULL,
  sl3_learner_OR = NULL,
  glm_formula_A = NULL,
  sl3_learner_A = NULL,
  glm_formula_Delta = NULL,
  sl3_learner_Delta = NULL,
  sl3_learner_default = Lrnr_hal9001_custom$new(max_degree = 2, smoothness_orders = 1,
    num_knots = c(20, 10))
)
}
\arguments{
\item{working_formula}{An working-model R formula object describing the functional form of the approximation of the conditional log odds ratio as a fnction of `W`.}

\item{W}{A named matrix of baseline covariates}

\item{A}{A binary vector with values in (0,1) encoding the treatment assignment}

\item{Y}{A binary outcome variable with values in (0,1)}

\item{Delta}{A binary vector that takes the value 1 if `Y` is osberved/not-missing and 0 otherwise.}

\item{weights}{An optional vector of weights for each observation. Use with caution. This can lead to invalid inferences if included naively.}

\item{W_new}{An optional matrix of new values of the baseline covariates `W` at which to predict odds ratio.}

\item{glm_formula_Y0W}{(Not recommended). An optional R formula object describing the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`}

\item{glm_formula_OR}{(Not recommended). An optional R formula object describing the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`}

\item{glm_formula_A}{(Not recommended). An optional R formula object describing the functional form of P(A=1|W). If provided, \code{glm} is used for the fitting. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)}

\item{sl3_learner_A}{An optional \code{tlverse/sl3} learner object used to estimate P(A=1|W).
If both \code{sl3_learner_A} and \code{glm_formula_A} are not provided, a default learner is used (Lrnr_hal9001).}

\item{glm_formula_Delta}{(Not recommended). An optional R formula object describing the functional form of P(Delta=1|A,W) to fit with glm. If provided, it is estimated using glm. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)}

\item{sl3_learner_Delta}{An optional \code{tlverse/sl3} learner object used to estimate P(Delta=1|A,W).}

\item{sl3_learner_default}{A default sl3 Learner to be used if neither a glm formula or sl3 learner is provided for one of the nuisance functions.
By default, Lrnr_hal9001 is used.}
}
\description{
Nonparametric Targeted estimates and inference for the conditional odds ratio with a user-specified parametric working-model with informative outcome missingness due to `A` and `W`.
This version also allows for the outcome to be missing-at-random conditional on A,W.
Note this version is totally nonparametric and requires no assumptions on the functional form of the conditional odds ratio for correct inference.
The user-specified working model is used to find a best-approximation (relative to the working-model) of the true conditional odds ratio.
This function can be viewed as a nonparametric version of the partially-linear logistic regression model where one uses the working model `logit(P(Y=1|A,W)) = b* A f(W) +  logit(P(Y=1|A=0,W)` for a user-specified parametric function `f(W)` to approximate the true function.
Nonparametrically correct and efficient inference is given for this best approximation of the conditional odds ratio.
}
