The Function
===============

```{coffee}
hpp <- function(data, lam = 0.3, phi = 0.2, bet1 = 0.1, delta = 0.6, 
                mu = 0, s = 1, a = 2, b = 2, v = 75, R = 50,
                bet.pr = matrix(c(0.1, 0.4, 0.1, 0.5), nrow = 4, ncol = 1),
                bet1.pr = matrix(c(0.2, 0.3), nrow = 2)){
  if((lam < 0) || (lam > 1))
    stop("lam must be between [0 and 1].")
  if((phi < 0) || (phi > 1))
    stop("phi must be between [0 and 1].")
  if((bet1 < 0) || (bet1 > 1))
    stop("bet1 must be between [0 and 1].")
  if((s < 0) || (s > 1))
    stop("s must be between [0 and 1].")
  if(a < 0)
    stop("a must be natural number.")
  if(b < 0)
    stop("b must be natural number.")
  if(!is.matrix(bet.pr))
    stop("bet.pr must be a 4 by 1 matrix.")
  if(!is.matrix(bet1.pr))
    stop("bet1.pr must be a 2 by 1 matrix.")
  
  n <- length(data)
  
  # The psy^{-1}
  psy <- 1 / (s ^ 2)
  
  # Compute the beta_0
  bet0 <- (1 - lam) * phi
  bet2 <- bet1 + delta
  
  # The beta parameter matrix
  bet <- matrix(c(bet0, bet1, bet0, bet2), nrow = 4, ncol = 1)
  
  PP <- numeric()
  
  x1 <- cbind(rep(1, v), data[1:v])
  x2 <- cbind(rep(1, n - v), data[(v + 1):n])
  o1 <- cbind(rep(0, v), rep(0, v))
  o2 <- cbind(rep(0, n - v), rep(0, n - v))
  
  # The X covariate matrix
  X <- rbind(cbind(x1, o1), cbind(o2, x2))

  # Case 1: Loop v from 1 to (n - 1)
  if(v < n){
    for(p in 1:R){
      
      # Generate random samples from N(mu, s) for innovation (errors)
      eps <- rnorm(n, mean = mu, sd = psy)
      eps1 <- c(0, eps[1:(n - 1)])
      
      # The u error matrix
      u <- cbind(eps - (lam * eps1))
      
      # Compute Z using the model below
      Z <- X %*% bet + u
      
      # Compute the sigma_{u}^{2}
      sigma.u <- psy * (1 + (lam ^ 2))
      
      # Compute the rho
      rho <- (-lam) / (1 + (lam ^ 2))
      
      Lambda <- matrix(NA, nrow = length(u), ncol = length(u))
      
      # Compute the diagonal matrix Lambda
      for(i in 1:length(u)){
        for(j in 1:length(u)){
          if(abs(i - j) == 0){
            Lambda[i, j] <- 1
          }
          if(abs(i - j) == 1){
            Lambda[i, j] <- rho
          }
          if(abs(i - j) > 1){
            Lambda[i, j] <- 0
          }
        }
      }
      
      # Compute the variance-covariance matrix
      v.cov <- sigma.u * Lambda
      
      # The identity matrix I_{4} (4 by 4)
      I <- matrix(0, nrow = ncol(X), ncol = ncol(X))
      diag(I) <- 1
      
      # The U matrix
      U <- t(X) %*% solve(Lambda) %*% X + I
      
      # The V matrix
      V <- t(X) %*% solve(Lambda) %*% Z + bet.pr
      
      # The W matrix
      W <- (2 * b) + t(Z) %*% solve(Lambda) %*% Z + (t(bet.pr) %*% bet.pr)
      
      # The M matrix
      M <- (-t(V)) %*% solve(U) %*% V + W
      
      # Compute the Posterior Probability
      PP[p] <- (det(Lambda)) ^ (- 1 / 2) * (2 / M) ^ (a + 2) * gamma(a + 2) * (det(U) ^ (1 / 2))
    }
  }
  if(v == n){
    for(p in 1:R){
      eps <- rnorm(n, mean = mu, sd = psy)
      eps1 <- c(0, eps[1:(n - 1)])
      
      # The u error matrix
      u <- cbind(eps - (lam * eps1))
      
      # Compute the sigma_{u}^{2}
      sigma.u <- psy * (1 + (lam ^ 2))
      
      # Compute the rho
      rho <- (-lam) / (1 + (lam ^ 2))
      
      Lambda <- matrix(NA, nrow = length(u), ncol = length(u))
      
      # Compute the diagonal matrix Lambda
      for(i in 1:length(u)){
        for(j in 1:length(u)){
          if(abs(i - j) == 0){
            Lambda[i, j] <- 1
          }
          if(abs(i - j) == 1){
            Lambda[i, j] <- rho
          }
          if(abs(i - j) > 1){
            Lambda[i, j] <- 0
          }
        }
      }
      
      # Compute the variance-covariance matrix
      v.cov <- sigma.u * Lambda
      
      # The B_{1} matrix
      B1 <- matrix(c(bet0, bet1), nrow = 2)
      
      # The X_{1} matrix
      X1 <- cbind(rep(1, n), data[1:n])
      
      # The Z_{1} matrix
      Z1 <- X1 %*% B1 + u
      
      # The identity matrix I_{2} (2 by 2)
      I2 <- matrix(0, nrow = 2, ncol = 2)
      diag(I2) <- 1
      
      # The U_{1} matrix
      U1 <- t(X1) %*% solve(Lambda) %*% X1 + I2
      
      # The V_{1} matrix
      V1 <- t(X1) %*% solve(Lambda) %*% Z1 + bet1.pr
      
      # The W_{1} matrix
      W1 <- (2 * b) + t(Z1) %*% solve(Lambda) %*% Z1 + (t(bet1.pr) %*% bet1.pr)
      
      # The M_{1} matrix
      M1 <- (-t(V1)) %*% solve(U1) %*% V1 + W1
      
      # Compute the new PP with v = n
      PP[p] <- (det(Lambda)) ^ (- 1 / 2) * (2 / M1) ^ (a + 1) * gamma(a + 1) * (det(U1) ^ (1 / 2))
    }
  }
  
  # Take the sum of all the posterior probability
  PP.sum <- sum(PP)
  
  # Take the proportion of the posterior probability
  PP.prop <- PP / PP.sum
  
  # Locate the point v at which Highest Posterior Probability (HPP) occurs
  max.v <- subset(1:R, PP.prop == max(PP.prop))
  
  # Return the output
  return(list("Posterior Probality" = PP.prop,
              "Highest Posterior Probability" = max(PP.prop),
              "The point v" = max.v))
}

hpp(data = airquality$Wind)
```
Misleads in this Function
===============
The function loops through the error term of the model Z = Xbeta + u, where u is the innovation (error). The function is supposed to estimate the location of the v in the time point of the series. But since the model is looped 50 times, thus the computation of the breakpoint is limited through the set of this domain only, that is 1 < v < 50. Hence if there are 100 time points in the data, then the computed HPP is absolutely misleading since v is computed up to 50 only, and 51 ahead are not included (therefore, BIAS).

Remedy
================
I did a remedy here, I made the looping through the error be general. And that is what the R = 50 argument meant in the function. So if the data has 100 samples, we can loop the process up to 100 as well by setting R = 100, to see if it can capture the true value of v.

Still BIAS
================
Though we did some remedy on the looping, but the parameter v, which suppose to be estimated is already set in the function. Which is of course bias. Should we not assign v to be unbias? But how will we assign the parameter X of Z = Xbeta + u. And if we assign any value v, we will have a new Z. What is the use of simulating a data, or using a realistic data, which suppose to be Z, if we there is new Z generated for every value of v?
