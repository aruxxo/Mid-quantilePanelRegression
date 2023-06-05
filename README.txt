This file describes the main "panelmidrq" function to perform mid-quantile regression for discrete panel responses.

It accompanies the paper "Mid-quantile regression for discrete panel data". 

## INPUT ARGUMENTS ##

formula : the specified model formula where "y" is the response variable; homogenous model "y ~ 1 + X"; fixed effects model "y ~ 1 + as.factor(ID) + X"; random effects model "y ~ 1 + (1|ID) + X"
data : the data.frame of interest with the response variable named "y"
TAU : a vector specifying the quantiles of interest 
cdf.model : a character indicating the model corresponding to the formula argument (either "HMG", "FE" or "RE")
lambda : a vector of pre-determined penalty terms 
lstar : logical; if TRUE the automatic penalty selection strategy described in the paper is performed to recover an optimal value for the shrinkage parameter
h : a numeric value for the link function; 0 is the identity link, 1 is the logarithmic link
stErr : logical; if TRUE the closed-form expression for standard errors calculation is applied and results are reported in the output

## STRUCTURE ## 
The function performs two distinct steps.
At the first step the conditional mid-CDF is estimated through a generalised linear mixed model approach where individual intercepts can be treated either as random parameters, as fixed coefficients or shrunk to a common intercept. At the second step, a possibly penalised objective function is maximised to obtain parameters estimates. A procedure based on the total variance law is implemented to calculate standard errors in closed form (when possible). The function includes also the procedure to select an automatic value for the penalty parameter.

## OUTPUT ##

Individual intercepts "alpha" and coefficients "beta" describing the effects of covariates are reported for al the quantiles-and-lambdas pairs of interest.
Standard errors are also quantile and penalty specific and reported in the array StErr (of corresponding dimensions). 
"Ghat" list collects objects involved in the estimation of the conditional mid-CDF at the first step.
"X" is the design model matrix.