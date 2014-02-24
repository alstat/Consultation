Determining the Breakpoint
============================

### Function for Simulating the Structural Change Data

```{coffee}
stchange.df <- function(data, mu = 0, s = sqrt(1), t = 20,
                        phi = 0.2, lam = 0.3, bet1 = 0.3, delta = 0.5){
  if (is.matrix(data)) {
    if (nrow(data) < 5) {
      stop('data must have sufficient time points.')
    } else {
      n <- nrow(data)
    }
  }
  if (is.numeric(data)) {
    if (length(data) < 5) {
      stop('data must have sufficient time points.')
    } else {
      n <- length(data)
    }
  }
  
  # Generate the error term
  error <- rnorm(n = n, mean = mu, sd = s)
  
  # The beta_2 parameter
  bet2 <- bet1 + delta
  
  # The beta_0 parameter
  bet0 <- (1 - lam) * phi
  
  # define the entries of Z matrix, with y1 from 1 to t
  # and y2 from n - t to t
  y1 <- matrix(0, nrow = t)
  y2 <- matrix(0, nrow = n - t)
  
  # Compute the entries of y1 and y2
  for(i in 2:t){ 
    y1[i, ] <- bet0 + bet1 * data[i - 1] + error[i] - 
      lam * (error[i - 1] - y1[(i - 1), ])
  }
  y2[1,] <- y1[t,]
  for(i in 2:(n - t)){
    y2[i, ] <- bet0 + bet2 * data[(t - 1) + (i - 1)] + error[t + i] - 
      lam * (error[t + i - 1] - y2[(i - 1), ])
  }
  
  # Combine the two variables to form the Z matrix
  Z <- rbind(y1, y2)
  
  # Return the simulated Z
  return('Structural Change Data' = Z)
}

```
### Function for Computing the Highest Posterior Probability

```{coffee}
hpp <- function(Z, Xcov, a = 2, b = 2, lam,
                bet.pr = matrix(c(0.1, 0.4, 0.1, 0.5), nrow = 4, ncol = 1),
                bet1.pr = matrix(c(0.2, 0.3), nrow = 2)){

  if (is.matrix(Z)) {
    if (nrow(Z) < 5) {
      stop('Z must have sufficient time points.')
    } else {
      n <- nrow(Z)
    }
  }
  if(is.numeric(Z)){
    if(length(Z) < 5){
      stop('Z must have sufficient time points.')
    } else {
      n <- length(Z)
    }
  }
  
  PP <- numeric()
  # Loop through v

  for(p in 2:n){
    x1 <- cbind(rep(1, p - 1), Xcov[1:(p - 1)])
    x2 <- cbind(rep(1, n - (p - 1)), Xcov[p:n])
    o1 <- cbind(rep(0, p - 1), rep(0, p - 1))
    o2 <- cbind(rep(0, n - (p - 1)), rep(0, n - (p - 1)))
    
    # The X covariate matrix
    X <- rbind(cbind(x1, o1), cbind(o2, x2))
    
    # Compute the rho
    rho <- (-lam) / (1 + (lam ^ 2))
    
    Lambda <- matrix(NA, nrow = n, ncol = n)
    
    # Compute the diagonal matrix Lambda
    for(i in 1:n){
      for(j in 1:n){
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
    # The identity matrix
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
  # The X_{1} covariate
  X1 <- cbind(rep(1, n), Xcov[1:n])
  
  # The Z_{1} matrix
  Z1 <- Z
  
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
  PP[p + 1] <- (det(Lambda)) ^ (- 1 / 2) * (2 / M1) ^ (a + 1) * gamma(a + 1) * (det(U1) ^ (1 / 2))

  # Obtain the complete cases of PP to remove NA
  PP <- PP[complete.cases(PP)]

  # Compute the sum of the Posterior Probability (PP)
  PP.sum <- sum(PP)

  # Take the proportion
  PP.prop <- PP/PP.sum  
  
  # Obtain the Highest Posterior Probability
  max.PP <- max(PP.prop)
  
  # Locate the Highest Posterior Probability
  v <- subset(1:n, PP.prop == max(PP.prop))
  
  # Return the output
  return(list('Posterior Probability' = PP.prop, 
              'Highest Posterior Probability' = max.PP, 
              'Location of v' = v))
  
}
````
### Example
Now let us take it into action, by using the **Wind** variable of the **airquality** data, with change occured at t = 100

```{coffee}
# Simulate structural change data
Z <- stchange.df(data = airquality$Wind, t = 100)

# Plot this data
plot.ts(Z)
```
![Plot](http://ubuntuone.com/6mK9NGjbb4yriELeZNsW0J)

Now, let us estimate the location of the change:
```{coffee}
hpp(Z, Xcov = airquality$Wind, lam = 0.3)

# OUTPUT
$`Posterior Probability`
  [1]           NA 1.076526e-05 1.937440e-05 3.472290e-05 4.478230e-05
  [6] 5.909355e-05 7.330088e-05 7.739827e-05 8.469058e-05 1.095821e-04
 [11] 1.121608e-04 1.210642e-04 1.292438e-04 1.379776e-04 1.452108e-04
 [16] 1.543793e-04 1.606308e-04 1.663876e-04 1.887898e-04 1.935662e-04
 [21] 1.992405e-04 2.054322e-04 2.171646e-04 2.208900e-04 2.273717e-04
 [26] 2.396598e-04 2.507611e-04 2.552600e-04 2.617741e-04 2.720536e-04
 [31] 2.815522e-04 2.982244e-04 3.116657e-04 3.263936e-04 3.440488e-04
 [36] 3.521945e-04 3.641926e-04 3.829103e-04 3.863730e-04 3.968397e-04
 [41] 4.076758e-04 4.208921e-04 4.341700e-04 4.510066e-04 4.666903e-04
 [46] 4.837754e-04 5.015541e-04 5.338420e-04 5.872963e-04 5.836855e-04
 [51] 6.011338e-04 6.149294e-04 6.305804e-04 6.672479e-04 7.159056e-04
 [56] 7.685259e-04 8.196898e-04 8.683809e-04 9.214306e-04 9.745537e-04
 [61] 1.034525e-03 1.037518e-03 1.062320e-03 1.127931e-03 1.166733e-03
 [66] 1.221159e-03 1.242114e-03 1.321236e-03 1.380764e-03 1.476736e-03
 [71] 1.551872e-03 1.687721e-03 1.851203e-03 2.043729e-03 2.274379e-03
 [76] 2.563755e-03 2.839054e-03 2.901513e-03 3.099275e-03 3.212228e-03
 [81] 3.445715e-03 3.983009e-03 4.207559e-03 4.583411e-03 5.244418e-03
 [86] 5.760919e-03 6.374002e-03 7.131844e-03 8.531848e-03 9.680368e-03
 [91] 1.079250e-02 1.241172e-02 1.441123e-02 1.703824e-02 2.422066e-02
 [96] 2.789608e-02 3.288710e-02 3.971390e-02 4.686181e-02 5.861113e-02
[101] 7.937039e-02 9.791611e-02 9.601788e-02 8.709169e-02 6.444179e-02
[106] 4.368289e-02 2.959643e-02 1.962886e-02 1.399850e-02 1.010635e-02
[111] 7.902020e-03 6.978552e-03 6.330590e-03 5.401910e-03 4.220773e-03
[116] 3.215448e-03 2.363690e-03 1.860073e-03 1.570190e-03 1.460716e-03
[121] 1.367811e-03 1.211395e-03 1.090836e-03 1.062215e-03 1.027700e-03
[126] 9.475844e-04 8.103569e-04 7.594271e-04 8.023569e-04 8.589711e-04
[131] 8.653655e-04 7.507137e-04 6.375399e-04 5.617624e-04 4.661041e-04
[136] 3.793584e-04 2.973394e-04 2.422137e-04 2.151591e-04 1.930329e-04
[141] 1.711122e-04 1.583191e-04 1.404009e-04 1.251847e-04 1.129757e-04
[146] 1.035608e-04 9.407195e-05 8.486119e-05 6.776356e-05 5.779964e-05
[151] 4.233256e-05 2.739350e-05 1.792025e-05 9.408889e-04

$`Highest Posterior Probability`
[1] 0.09791611

$`Location of v`
[1] 102
```
Since Z was generated using lam = 0.3, and we estimated the v using the same value of lam. So to make it not bias (unless it was coincidental), lets assign it to lam = 0.5, here is the estimated v.
```{coffee}
hpp(Z, Xcov = airquality$Wind, lam = 0.5)

# OUTPUT
$`Posterior Probability`
  [1]           NA 2.662452e-06 5.209650e-06 9.923727e-06 1.356344e-05
  [6] 1.861544e-05 2.382400e-05 2.528947e-05 2.806716e-05 3.810424e-05
 [11] 3.859087e-05 4.236729e-05 4.602576e-05 5.026326e-05 5.320055e-05
 [16] 5.530818e-05 5.780167e-05 6.023789e-05 6.837555e-05 7.046195e-05
 [21] 7.287500e-05 7.665408e-05 8.163130e-05 8.443892e-05 8.783298e-05
 [26] 9.455981e-05 1.006575e-04 1.023722e-04 1.057036e-04 1.113456e-04
 [31] 1.165890e-04 1.264974e-04 1.367184e-04 1.464654e-04 1.542042e-04
 [36] 1.583649e-04 1.654715e-04 1.725506e-04 1.773314e-04 1.865829e-04
 [41] 1.931327e-04 1.991297e-04 2.053736e-04 2.119605e-04 2.231212e-04
 [46] 2.363856e-04 2.475428e-04 2.642187e-04 3.116860e-04 3.114255e-04
 [51] 3.213134e-04 3.295630e-04 3.358664e-04 3.607173e-04 3.903841e-04
 [56] 4.226611e-04 4.555880e-04 4.858791e-04 5.190691e-04 5.561718e-04
 [61] 6.042616e-04 6.111246e-04 6.243874e-04 6.570276e-04 6.957887e-04
 [66] 7.389677e-04 7.556906e-04 8.127534e-04 8.378476e-04 8.839311e-04
 [71] 9.461078e-04 1.027895e-03 1.129388e-03 1.299217e-03 1.572467e-03
 [76] 1.950137e-03 2.292309e-03 2.349699e-03 2.519180e-03 2.623055e-03
 [81] 2.776754e-03 3.320818e-03 3.644154e-03 4.231022e-03 5.118490e-03
 [86] 5.632186e-03 6.167580e-03 7.047759e-03 8.843860e-03 9.856748e-03
 [91] 1.123828e-02 1.289109e-02 1.590700e-02 1.897586e-02 2.594479e-02
 [96] 2.879213e-02 3.225415e-02 3.751111e-02 4.271119e-02 5.285125e-02
[101] 8.336169e-02 1.149278e-01 1.154094e-01 9.875210e-02 7.305029e-02
[106] 4.451027e-02 2.811176e-02 1.734446e-02 1.122242e-02 7.983040e-03
[111] 6.113836e-03 4.921205e-03 4.240502e-03 3.319909e-03 2.483997e-03
[116] 1.763566e-03 1.286351e-03 1.025022e-03 8.628734e-04 7.829694e-04
[121] 7.041350e-04 6.360812e-04 5.880786e-04 5.549800e-04 5.193727e-04
[126] 4.606261e-04 3.847560e-04 3.609334e-04 3.985204e-04 4.144796e-04
[131] 4.046650e-04 3.538237e-04 3.008334e-04 2.650294e-04 2.107970e-04
[136] 1.739831e-04 1.371163e-04 1.107368e-04 9.500630e-05 8.416625e-05
[141] 7.303587e-05 6.682833e-05 5.993010e-05 5.291312e-05 4.691950e-05
[146] 4.089457e-05 3.511456e-05 3.060220e-05 2.308396e-05 1.918119e-05
[151] 1.311278e-05 7.383333e-06 4.390280e-06 4.996218e-04

$`Highest Posterior Probability`
[1] 0.1154094

$`Location of v`
[1] 103
```
So we have a good estimate, and base on the plot, it is in 102 where the big change occur.

### HPP R = 50 runs

```{coffee}
hpp.runs <- function(dat, mu = 0, s = sqrt(1), t = 20, phi = 0.2, lam = 0.3, 
                     bet1 = 0.3, delta = 0.5, a = 2, b = 2, R = 50, 
                     bet.pr = matrix(c(0.1, 0.4, 0.1, 0.5), nrow = 4, ncol = 1),
                     bet1.pr = matrix(c(0.2, 0.3), nrow = 2)){
  hpp.out <- matrix(NA, nrow = R, ncol = 2)
  colnames(hpp.out) <- c('HPP at v', 'HPP near v')
  
  # Loop the process R times
  for(i in 1:R){
    # Generate the data
    d <- stchange.df(dat, mu = mu, s = s, t = t, phi = phi, 
                     lam = lam, bet1 = bet1, delta = delta)
    
    # Compute the HPP
    h <- hpp(d, Xcov = dat, a = a, b = b, lam = lam,
             bet.pr = bet.pr, bet1.pr = bet1.pr)
    
    # Extract the location of v
    h.out <- h$'Location of v'
    
    # Compute the HPP - 5% of R runs
    h.nvl <- h.out - (0.05 * R)
    
    # Compute the HPP + 5% of R runs
    h.nvu <- h.out + (0.05 * R)
    
    # Obtain the HPP at v
    if(h.out == t){
      hpp.out[i, 1] <- 1
    }
    
    # Obtain the HPP near v
    if((t > h.nvl) && (t < h.nvu)){
      hpp.out[i, 2] <- 1
    }
  }
  
  # Return the output
  return(list('HPP at v' = sum(hpp.out[,1], na.rm = TRUE),
              'HPP near v' = sum(hpp.out[,2], na.rm = TRUE),
              'Percentage near v' = sum(hpp.out[,2], na.rm = TRUE)/R,
              'HPP Runs' = hpp.out))
}
```
Let's give it a try,
```{coffee}
hpp.runs(dat = airquality$Wind[1:20], lam = 0.3, t = 8, R = 50)

# OUTPUT
$`HPP at v`
[1] 0

$`HPP near v`
[1] 50

$`Percentage near v`
[1] 1

$`HPP Runs`
      HPP at v HPP near v
 [1,]       NA          1
 [2,]       NA          1
 [3,]       NA          1
 [4,]       NA          1
 [5,]       NA          1
 [6,]       NA          1
 [7,]       NA          1
 [8,]       NA          1
 [9,]       NA          1
[10,]       NA          1
[11,]       NA          1
[12,]       NA          1
[13,]       NA          1
[14,]       NA          1
[15,]       NA          1
[16,]       NA          1
[17,]       NA          1
[18,]       NA          1
[19,]       NA          1
[20,]       NA          1
[21,]       NA          1
[22,]       NA          1
[23,]       NA          1
[24,]       NA          1
[25,]       NA          1
[26,]       NA          1
[27,]       NA          1
[28,]       NA          1
[29,]       NA          1
[30,]       NA          1
[31,]       NA          1
[32,]       NA          1
[33,]       NA          1
[34,]       NA          1
[35,]       NA          1
[36,]       NA          1
[37,]       NA          1
[38,]       NA          1
[39,]       NA          1
[40,]       NA          1
[41,]       NA          1
[42,]       NA          1
[43,]       NA          1
[44,]       NA          1
[45,]       NA          1
[46,]       NA          1
[47,]       NA          1
[48,]       NA          1
[49,]       NA          1
[50,]       NA          1
```

