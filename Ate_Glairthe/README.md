## Error Growth Algorithm
Below is the function for computing the statistical inference of unperturbed, perturbed and the error. The parameter of the function, `x`, is the series/data points generated from the models.
```{coffee}
stat_moments <- function(x) {
  if (!require(moments)) 
    install.packages("moments")
  library(moments)
  
  if (is.matrix(x)) {
    output <- matrix(NA, nrow = 5, ncol = ncol(x))
    colnames(output) <- c("Unperturbed", "Perturbed")
    rownames(output) <- c("Mean", "Std Dev", "Variance", "Kurtosis", "Skewness")
    
    for (i in 1:ncol(x)) {
      output[, i] <- c(mean(x[, i]), sd(x[, i]), var(x[, i]), 
                    kurtosis(x[, i]), skewness(x[, i]))
    }  
  } else if (is.numeric(x)) {
    output <- c(mean(x), sd(x), var(x), 
             kurtosis(x), skewness(x))
    names(output) <- c("Mean", "Std Dev", "Variance", "Kurtosis", "Skewness")
  }
    
  return(as.data.frame(output))
}
```
Another function below computes the error growth of the unperturbed and perturbed model. The parameter of the function which are the following:

1.  `n` - the number of initial points generated, default is set to 10;
2.  `r` - the number of iterations needed for simulating the z series;
3.  `min` - the minimum parameter of the uniform distribution for generating the initial points;
4.  `max` - the maximum parameter of the uniform distribution for generating the initial points;
5.  `alpha` - the alpha parameter, set to 1.33;
6.  `epsilon` - the level of magnitude, default to 10 ^ (-9);
7.  `yt_min` - the minimum parameter of the uniform distribution for generating the errors `yt`;
8.  `yt_max` - the maximum parameter of the uniform distribution for generating the errors `yt`;
9.  `output` - if set to `"partial"`, then only first 10 of the iteration outputs are shown, especially for z series generated through `r` iterations. Else if set to `"full"`, then all output are returned.

```{coffee}
error_growth <- function (n = 10, r = 100, min = -2, max = 2, alpha = 1.33, 
                          epsilon = 10 ^ (-9), yt_min = 0, yt_max = 1, 
                          output = "partial") {
  # Generate n (x, y, z) from uniform distribution
  # with a = -2 and b = 2:
  x <- x1 <- runif(n, min = min, max = max)
  y <- y1 <- runif(n, min = min, max = max)
  z <- z1 <- runif(n, min = min, max = max)
  
  # Initialize
  z_vec <- matrix(NA, nrow = r, ncol = 2)
  colnames(z_vec) <- c("Unperturbed", "Perturbed")
  z_vec_list <- z_vec_moments <- z_errors_list <- z_errors_moments <- list()
  yt <- runif(n, min = yt_min, max = yt_max)
  
  for (i in 1:n) {
    # For every initial point, simulate the unperturbed of 
    # the Henon type model with 100,000 iterates.
    for (j in 1:r) {
      x[i] <- 1.4 - (x[i] ^ 2) + (0.3 * y[i])
      y[i] <- x[i]
      z[i] <- (alpha * x[i] * sin(2 * pi * z[i])) / 2 * pi
      z_vec[j, 1] <- z[i]
    }
    x[i + 1] <- x[i]
    y[i + 1] <- x[i]
    z[i + 1] <- z[i]
    
    # Repeat step 2 and 3, but consider the perturbed model of 
    # the Henon type model with error level magnitude at 10 ^ (-9)
    # from uniform distribution (0, 1).
    for (j in 1:r) {
      x1[i] <- 1.4 - (x1[i] ^ 2) + (0.3 * y1[i])
      y1[i] <- x1[i]
      z1[i] <- (alpha * x1[i] * sin(2 * pi * z1[i])) / 2 * pi
      z_vec[j, 2] <- z1[i] + (epsilon * yt[i])
    }
    x1[i + 1] <- x1[i]
    y1[i + 1] <- x1[i]
    z1[i + 1] <- z1[i]
    
    # Extract the list of the simulated zs
    z_vec_list[[i]] <- as.data.frame(z_vec)
    
    # Get the statistical averages (mean, 
    # variance and other moments) of z of the simulated model.
    z_vec_moments[[i]] <- stat_moments(z_vec)
    
    # Get the error of the perturbed and unperturbed model.
    z_errors <- z_vec_list[[i]][, 1] - z_vec_list[[i]][, 2]
    z_errors_list[[i]] <- as.data.frame(z_errors)
    z_errors_moments[[i]] <- stat_moments(z_errors)
  }
  
  if (output == "partial") {
    z_vec_list <- lapply(z_vec_list, function(x) head(x))
    z_errors_list <- lapply(z_errors_list, function(x) head(x))
  } else if (output == "full") {
    z_vec_list <- z_vec_list; z_errors_list <- z_errors_list
  }
  
  return(list(z_list = z_vec_list, z_moments = z_vec_moments, z_errors = z_errors_list,
         z_errors_moments = z_errors_moments))
}
```
Lets try to run it, using default values we have
```{coffee}
error_growth()

## OUTPUT

$z_list
$z_list[[1]]
  Unperturbed   Perturbed
1 -2.00343286 -2.00343285
2 -0.01942366 -0.01942366
3 -0.34169581 -0.34169581
4  0.00337743  0.00337743
5  0.06203725  0.06203725
6 -0.10998701 -0.10998701

$z_list[[2]]
  Unperturbed  Perturbed
1   2.2357056  2.2357056
2  -0.3740166 -0.3740166
3  -1.9527229 -1.9527229
4   0.0416762  0.0416762
5   0.7657149  0.7657149
6   0.3737030  0.3737030

$z_list[[3]]
  Unperturbed   Perturbed
1 -1.76305052 -1.76305052
2 -0.37427365 -0.37427365
3 -1.94960324 -1.94960323
4  0.04433921  0.04433921
5  0.81340971  0.81340971
6  0.34612243  0.34612243

$z_list[[4]]
  Unperturbed  Perturbed
1   0.6917247  0.6917247
2   0.3506416  0.3506416
3   2.2139529  2.2139529
4   0.1387560  0.1387560
5   2.2642450  2.2642450
6  -0.3740323 -0.3740323

$z_list[[5]]
  Unperturbed  Perturbed
1  -1.8856735 -1.8856735
2  -0.2471518 -0.2471518
3  -2.7442185 -2.7442185
4   0.1422987  0.1422987
5   2.3060435  2.3060435
6  -0.3524924 -0.3524924

$z_list[[6]]
  Unperturbed   Perturbed
1  0.73878759  0.73878759
2  0.37460393  0.37460393
3  1.94559004  1.94559004
4 -0.04773676 -0.04773676
5 -0.87392758 -0.87392758
6 -0.26732692 -0.26732692

$z_list[[7]]
  Unperturbed   Perturbed
1  1.74295306  1.74295306
2  0.37516741  0.37516741
3  1.93872378  1.93872378
4 -0.05347815 -0.05347815
5 -0.97527462 -0.97527462
6 -0.05810660 -0.05810660

$z_list[[8]]
  Unperturbed  Perturbed
1  -2.1902139 -2.1902139
2   0.3493495  0.3493495
3   2.2270498  2.2270498
4   0.1409148  0.1409148
5   2.2898501  2.2898501
6  -0.3638251 -0.3638251

$z_list[[9]]
  Unperturbed   Perturbed
1  0.72799633  0.72799633
2  0.37195220  0.37195220
3  1.97757365  1.97757365
4 -0.01999805 -0.01999805
5 -0.37067920 -0.37067920
6  0.27265402  0.27265402

$z_list[[10]]
  Unperturbed  Perturbed
1   2.3025365  2.3025365
2  -0.3552607 -0.3552607
3  -2.1659467 -2.1659467
4  -0.1229923 -0.1229923
5  -2.0649562 -2.0649562
6   0.1490482  0.1490482


$z_moments
$z_moments[[1]]
         Unperturbed   Perturbed
Mean     -0.08381530 -0.08381530
Std Dev   1.40684384  1.40684384
Variance  1.97920960  1.97920960
Kurtosis  2.70546754  2.70546754
Skewness  0.00686246  0.00686246

$z_moments[[2]]
         Unperturbed  Perturbed
Mean      -0.1760476 -0.1760476
Std Dev    1.3082051  1.3082051
Variance   1.7114005  1.7114005
Kurtosis   2.8490025  2.8490025
Skewness   0.1695342  0.1695342

$z_moments[[3]]
         Unperturbed  Perturbed
Mean      -0.1803717 -0.1803717
Std Dev    1.4003994  1.4003994
Variance   1.9611186  1.9611186
Kurtosis   2.5662690  2.5662690
Skewness   0.2145500  0.2145500

$z_moments[[4]]
         Unperturbed  Perturbed
Mean       0.1260126  0.1260126
Std Dev    1.3980872  1.3980872
Variance   1.9546477  1.9546477
Kurtosis   2.5304314  2.5304314
Skewness  -0.1441903 -0.1441903

$z_moments[[5]]
         Unperturbed   Perturbed
Mean      0.06777242  0.06777242
Std Dev   1.45753840  1.45753840
Variance  2.12441818  2.12441818
Kurtosis  2.38848930  2.38848930
Skewness -0.10567715 -0.10567715

$z_moments[[6]]
         Unperturbed  Perturbed
Mean      -0.2805551 -0.2805551
Std Dev    1.3250029  1.3250029
Variance   1.7556327  1.7556327
Kurtosis   2.7009547  2.7009547
Skewness   0.2944671  0.2944671

$z_moments[[7]]
         Unperturbed   Perturbed
Mean     -0.02960195 -0.02960195
Std Dev   1.37642574  1.37642574
Variance  1.89454781  1.89454781
Kurtosis  2.60695680  2.60695680
Skewness -0.01464343 -0.01464343

$z_moments[[8]]
         Unperturbed   Perturbed
Mean      0.01489319  0.01489319
Std Dev   1.39433977  1.39433977
Variance  1.94418339  1.94418339
Kurtosis  2.59947938  2.59947938
Skewness -0.13879751 -0.13879751

$z_moments[[9]]
         Unperturbed  Perturbed
Mean      0.10029736 0.10029736
Std Dev   1.34876088 1.34876088
Variance  1.81915590 1.81915590
Kurtosis  2.63116432 2.63116432
Skewness  0.00673927 0.00673927

$z_moments[[10]]
         Unperturbed   Perturbed
Mean     -0.01451972 -0.01451972
Std Dev   1.34366274  1.34366274
Variance  1.80542955  1.80542955
Kurtosis  2.65327230  2.65327230
Skewness  0.08645403  0.08645403


$z_errors
$z_errors[[1]]
       z_errors
1 -5.335763e-10
2 -5.335764e-10
3 -5.335764e-10
4 -5.335764e-10
5 -5.335764e-10
6 -5.335764e-10

$z_errors[[2]]
       z_errors
1 -6.882681e-10
2 -6.882682e-10
3 -6.882681e-10
4 -6.882682e-10
5 -6.882682e-10
6 -6.882682e-10

$z_errors[[3]]
      z_errors
1 -6.96359e-10
2 -6.96359e-10
3 -6.96359e-10
4 -6.96359e-10
5 -6.96359e-10
6 -6.96359e-10

$z_errors[[4]]
       z_errors
1 -8.030511e-10
2 -8.030511e-10
3 -8.030510e-10
4 -8.030511e-10
5 -8.030510e-10
6 -8.030511e-10

$z_errors[[5]]
       z_errors
1 -8.986767e-10
2 -8.986767e-10
3 -8.986767e-10
4 -8.986767e-10
5 -8.986767e-10
6 -8.986766e-10

$z_errors[[6]]
       z_errors
1 -5.284952e-10
2 -5.284952e-10
3 -5.284952e-10
4 -5.284952e-10
5 -5.284952e-10
6 -5.284952e-10

$z_errors[[7]]
       z_errors
1 -2.420248e-10
2 -2.420249e-10
3 -2.420248e-10
4 -2.420249e-10
5 -2.420250e-10
6 -2.420249e-10

$z_errors[[8]]
       z_errors
1 -2.883476e-10
2 -2.883476e-10
3 -2.883476e-10
4 -2.883475e-10
5 -2.883476e-10
6 -2.883476e-10

$z_errors[[9]]
       z_errors
1 -4.109650e-10
2 -4.109651e-10
3 -4.109650e-10
4 -4.109651e-10
5 -4.109651e-10
6 -4.109651e-10

$z_errors[[10]]
       z_errors
1 -3.499725e-10
2 -3.499725e-10
3 -3.499725e-10
4 -3.499725e-10
5 -3.499725e-10
6 -3.499725e-10


$z_errors_moments
$z_errors_moments[[1]]
                output
Mean     -5.335764e-10
Std Dev   6.913266e-17
Variance  4.779324e-33
Kurtosis  2.176964e+00
Skewness  2.321293e-01

$z_errors_moments[[2]]
                output
Mean     -6.882682e-10
Std Dev   4.947035e-17
Variance  2.447316e-33
Kurtosis  1.210200e+00
Skewness  3.229036e-01

$z_errors_moments[[3]]
                output
Mean     -6.963590e-10
Std Dev   5.169558e-18
Variance  2.672433e-35
Kurtosis  3.740651e+00
Skewness -1.633647e+00

$z_errors_moments[[4]]
                output
Mean     -8.030510e-10
Std Dev   5.329679e-17
Variance  2.840547e-33
Kurtosis  1.143607e+00
Skewness  3.474917e-01

$z_errors_moments[[5]]
                output
Mean     -8.986767e-10
Std Dev   2.229346e-17
Variance  4.969984e-34
Kurtosis  1.433784e+00
Skewness  2.750495e-01

$z_errors_moments[[6]]
                output
Mean     -5.284952e-10
Std Dev   9.086157e-17
Variance  8.255825e-33
Kurtosis  3.027597e+00
Skewness  1.423878e+00

$z_errors_moments[[7]]
                output
Mean     -2.420249e-10
Std Dev   7.794859e-17
Variance  6.075982e-33
Kurtosis  2.225239e+00
Skewness -8.157583e-01

$z_errors_moments[[8]]
                output
Mean     -2.883476e-10
Std Dev   1.254596e-17
Variance  1.574010e-34
Kurtosis  1.505620e+00
Skewness  6.734441e-01

$z_errors_moments[[9]]
                output
Mean     -4.109651e-10
Std Dev   2.591459e-17
Variance  6.715659e-34
Kurtosis  1.104623e+00
Skewness  6.933973e-02

$z_errors_moments[[10]]
                output
Mean     -3.499725e-10
Std Dev   1.051669e-17
Variance  1.106008e-34
Kurtosis  2.362867e+00
Skewness  9.600309e-01
```
