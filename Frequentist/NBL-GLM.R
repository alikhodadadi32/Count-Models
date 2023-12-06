# This is a high-level documentation, and specific details of the 
# Maximum Likelihood Calculations and the negative binomial Lindley
# distribution implementation may need to
# be adapted based on your specific needs.
# Consider adding more information about the distribution itself,
# its parameters, the initai values for the maximization algorithms, 
# and the intended use of the MLE estimates.


# This function calculates the negative binomial Lindley (NBL) likelihood 
# function for a set of covariates and crash data.

# Inputs:
# v: A vector of length 4 containing the following parameters:
# beta_0: The intercept in the log-linear model.
# beta_1: The coefficient for the log of AADT in the log-linear model.
# beta_2: The coefficient for the road length in the log-linear model.
# r: The overdispersion parameter for the negative binomial distribution.

# data: A data frame containing the following columns:
#     AADT: The average annual daily traffic (AADT) for each road segment.
#     Length: The length of each road segment.
#     Crash: The number of crashes observed on each road segment.

# Outputs:
# ll: The negative log-likelihood of the observed crash data given the model and parameter values.


MLE <- function(data, par){
  beta_0 <- par[1]
  beta_1 <- par[2]
  beta_2 <- par[3]
  r <- par[4]

  ll <- 0
  for (i in 1:length(data["Crash"])){
    A <- 0
    theta <- exp(beta_0 + beta_1 * data["AADT"][i] + beta_2 * data["Length"][i])
    for (l in 0:data["Crash"][i]){
      A <- A + choose(data[i],l) * ((-1)^l) * 
      ((theta + r + l + 1)/(theta + r + l)^2)
    }
    ll <- ll + 2*log(theta)-log(theta+1)+log((A)) + lfactorial(r+data["Crash"][i]-1)-lfactorial(r-1)-lfactorial(data["Crash"][i])   
  }
  return(-ll)
}


# solution (1) for maximizing the likelihood funciton: 

# This part defines a genetic algorithm (GA) to estimate the maximum likelihood (MLE) of the negative binomial Lindley distribution parameters.
# Parameters:
# type: Specifies the type of parameters to be optimized. Set to "real-valued" for continuous parameters.
# fitness: Defines the fitness function. Here, it is set to MLE.
# lower: Lower bounds for the parameters.
# upper: Upper bounds for the parameters.
# maxiter: Maximum number of iterations.
# parallel: Run the GA in parallel mode (Boolean).
# popSize: Population size for the GA.

# This code assumes the availability of a function MLE that calculates the MLE for the negative binomial Lindley distribution given a vector of data and covariates (v).
# The number of distribution parameters estimated must match the length of the lower and upper vectors.

GA <- ga(type = "real-valued", fitness = MLE , lower = c(-100,-100,-100,-100),
upper = c(100,100,100,100), maxiter = 500 , parallel = F, popSize = 200)


# solution (2) for maximizing the likelihood funciton:

# optim(par = (0,0,0,0) , fn = MLE ,method = "L-BFGS-B",
#  lower = c(-100,-100,-100,-100), upper = c(100,100,100,100))
