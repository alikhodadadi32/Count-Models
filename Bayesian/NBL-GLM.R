# Negative Binomial Lindley Generalized Linear Model (NBL-GLM) with JAGS
# refer to -- https://www.sciencedirect.com/science/article/pii/S0001457521001342 -- for more information about the model  

# Introduction:
# The Negative Binomial Lindley Generalized Linear Model (NBL-GLM) is a statistical model used to analyze overdispersed count data, often encountered in ecological and epidemiological studies.

# Dependencies:
# The R code provided in this repository requires either one of the following R packages:
library("R2WinBUGS")
library("R2jags")



# Model Specification:
# The NBL-GLM model is specified using a text file named "NBL-GLM.txt". The model specification includes the following components:
# Likelihood: The likelihood function for the Negative Binomial distribution with Lindley-adjusted mean function.
# Priors: Prior distributions are specified for the model parameters, including the regression coefficients (beta_0, beta_1, beta_2, beta_3), the overdispersion parameter (phi), and the Lindley parameter (theta).
# The model assumes a Negative Binomial Lindley distribution for the response variable y, with additional components to model overdispersion, frality term, and incorporate covariates.
# The sink function is used to redirect the model specification to a text file.
sink ("NBL-GLM.txt")
  cat("
  model {
      for(i in 1:N) {
        y[i] ~ dnegbin(p[i],phi)                              
        p[i] <- (phi)/((phi) + adj_mu[i])
        mu[i] <- exp(beta_0 + beta_1 * x_1[i] + beta_2 * x_2[i] + beta_3 * x_3[i])
        adj_mu[i] <- mu[i] * eps[i]
        z[i]~ dbern(Theta)
        eps[i] ~ dgamma(1 + z [i], theta)
        ll[i] <- logdensity.negbin(y[i], p[i], phi)
      }
    beta_0 ~ dnorm(0,0.001)
    beta_1 ~ dnorm(0,0.001)
    beta_2 ~ dnorm(0,0.001)
    beta_3 ~ dnorm(0,0.001)
    Theta ~ dbeta(N/3,N/2)
    theta <-  (1/Theta)-1
    phi ~ dgamma(0.1,0.1)
  }
  ", fill = TRUE)
  sink()


# The MCMC process is run using the jags function. The jags function provides various options for controlling the MCMC process. 
# The following parameters are specified:
# data: A list containing the data for the model.
# parameters.to.save: A vector of parameter names to save from the MCMC chains.
# model.file: The filename of the JAGS model specification file.
# n.chains: The number of MCMC chains to run.
# n.iter: The total number of MCMC iterations.
# n.burnin: The number of burn-in iterations.
# n.thin: The thinning factor for saving MCMC samples.
jags(data, parameters.to.save, model.file="NB-GLM.txt", n.chains, n.iter, n.burnin, n.thin)
