Generating Structural Change Time Series
============================

The model use for generating the structural data is,

Z = X\beta + u

In R, this can be coded into the following

```{coffee}
stchange.df <- function(data, lam = 0.3, phi = 0.2, bet1 = 0.1, 
                        delta = 0.5, mu = 0, s = 1, v = 3){
  
  if((lam < 0) || (lam > 1))
    stop("lam must be between [0 and 1].")
  if((phi < 0) || (phi > 1))
    stop("phi must be between [0 and 1].")
  if((bet1 < 0) || (bet1 > 1))
    stop("bet1 must be between [0 and 1].")
  if((s < 0) || (s > 1))
    stop("s must be between [0 and 1].")
  
  # The psy^{-1}
  psy <- 1 / (s ^ 2)
  
  # Compute the beta_0
  bet0 <- (1 - lam) * phi
  
  # The beta parameter matrix
  bet2 <- bet1 + delta
  bet <- matrix(c(bet0, bet1, bet0, bet2), nrow = 4, ncol = 1)
  
  # Generate random samples from N(mu, s) for innovation (errors)
  n <- length(data)
  eps <- rnorm(n, mean = mu, sd = psy)
  eps1 <- c(0, eps[1:(n - 1)])
  
  # The u error matrix
  u <- cbind(eps - (lam * eps1))
  
  PP <- numeric()
  
  # Case 1: Loop v from 1 to (n - 1)
  x1 <- cbind(rep(1, v), data[1:v])
  x2 <- cbind(rep(1, n - v), data[(v + 1):n])
  o1 <- cbind(rep(0, v), rep(0, v))
  o2 <- cbind(rep(0, n - v), rep(0, n - v))
  
  # The X covariate matrix
  X <- rbind(cbind(x1, o1), cbind(o2, x2))
  
  # Compute Z using the model below
  Z <- X %*% bet + u
  
  return(Z)
}
```
So to generate the data using **Wind** variable of **airquality**, we simply run.
```{coffee}
a <- stchange.df(data = airquality$Wind, v = 75, delta = 0.6)
```
Plotting this, gives us
```{coffee}
plot.ts(a)
```
![plot](http://ubuntuone.com/1f7y9o69CzvoosJ7ZkCpbr)

It is clear that there is a structural change in the data, but at what point did the change occur?
