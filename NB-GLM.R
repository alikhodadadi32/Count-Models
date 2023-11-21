library("R2WinBUGS")
library("R2jags")


# The code below makes the NB.txt file which contains the Negative Binomial Generalized Linear Model (NB-GLM) model specification required by JAGS. 
# y: dependent variable or counts, X_: covariates, N: number of total observation
# refer to -- https://www.sciencedirect.com/science/article/pii/S0001457521001342 -- for more information about the model  

sink ("NB-GLM.txt")
  cat("
  model {
      for(i in 1:N) {
        y[i] ~ dnegbin(p[i],phi)                              
        p[i] <- (phi)/((phi)+mu[i])
        mu[i] <- exp(beta_0 + beta_1 * x_1[i] + beta_2 * x_2[i] + beta_3 * x_3[i])
        ll[i] <- logdensity.negbin(y[i], p[i], phi)
      }
    beta_0 ~ dnorm(0,0.001)
    beta_1 ~ dnorm(0,0.001)
    beta_2 ~ dnorm(0,0.001)
    beta_3 ~ dnorm(0,0.001)
    phi ~ dgamma(0.1,0.1)
  }
  ", fill = TRUE)
  sink()


# Run the MCMC process
jags(data, parameters.to.save, model.file="NB.txt", n.chains, n.iter, n.burnin, n.thin)
