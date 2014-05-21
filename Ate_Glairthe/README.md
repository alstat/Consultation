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
  
  return(list("z series from iteration" = z_vec_list, 
              "Statistical inference of z series from iteration" = z_vec_moments, 
              "errors of z series from iteration" = z_errors_list,
              "Statistical inference of errors from iteration" = z_errors_moments))
}
```
Lets try to run it, using default values we have
```{coffee}
error_growth()

## OUTPUT

$`z series from iteration`
$`z series from iteration`[[1]]
  Unperturbed   Perturbed
1 -3.41517722 -3.41517722
2  0.98871744  0.98871744
3 -0.03742858 -0.03742858
4 -0.68730508 -0.68730508
5 -0.32776829 -0.32776829
6 -2.43527207 -2.43527207

$`z series from iteration`[[2]]
  Unperturbed   Perturbed
1 -0.37141035 -0.37141035
2 -1.98404204 -1.98404204
3  0.01425327  0.01425327
4  0.26453770  0.26453770
5 -0.37396958 -0.37396958
6 -1.95329074 -1.95329074

$`z series from iteration`[[3]]
  Unperturbed  Perturbed
1 -0.17957411 -0.1795741
2 -2.48030474 -2.4803047
3 -0.01757601 -0.0175760
4 -0.32597999 -0.3259800
5  0.33354843  0.3335484
6  2.37508668  2.3750867

$`z series from iteration`[[4]]
  Unperturbed  Perturbed
1  -0.2513374 -0.2513374
2  -2.7445611 -2.7445611
3   0.1423095  0.1423095
4   2.3061693  2.3061693
5  -0.3523900 -0.3523900
6  -2.1959992 -2.1959992

$`z series from iteration`[[5]]
  Unperturbed  Perturbed
1   0.1139636  0.1139636
2   1.8016299  1.8016299
3  -0.1349657 -0.1349657
4  -2.2182838 -2.2182838
5   0.3681035  0.3681035
6   2.0230156  2.0230156

$`z series from iteration`[[6]]
   Unperturbed    Perturbed
1  0.369919335  0.369919335
2  2.001721571  2.001721572
3  0.001540224  0.001540224
4  0.028624038  0.028624039
5 -0.067176582 -0.067176582
6 -1.124380579 -1.124380578

$`z series from iteration`[[7]]
   Unperturbed    Perturbed
1 -0.002610316 -0.002610315
2 -0.045013392 -0.045013391
3 -0.039737781 -0.039737780
4 -0.730862609 -0.730862608
5 -0.372823879 -0.372823878
6 -1.967120073 -1.967120072

$`z series from iteration`[[8]]
  Unperturbed  Perturbed
1  -0.2868242 -0.2868242
2  -2.6715193 -2.6715193
3   0.1254288  0.1254288
4   2.0971323  2.0971323
5  -0.2152244 -0.2152244
6  -2.6793994 -2.6793994

$`z series from iteration`[[9]]
  Unperturbed   Perturbed
1 -0.37541227 -0.37541227
2 -1.93573248 -1.93573248
3  0.05594890  0.05594890
4  1.01850545  1.01850545
5 -0.04356637 -0.04356637
6 -0.74196283 -0.74196283

$`z series from iteration`[[10]]
  Unperturbed   Perturbed
1 -0.32518675 -0.32518675
2 -2.44404437 -2.44404437
3 -0.04903734 -0.04903734
4 -0.89698982 -0.89698982
5 -0.22644058 -0.22644058
6 -2.71464191 -2.71464191


$`Statistical inference of z series from iteration`
$`Statistical inference of z series from iteration`[[1]]
         Unperturbed  Perturbed
Mean      -0.1763497 -0.1763497
Std Dev    1.3492610  1.3492610
Variance   1.8205051  1.8205051
Kurtosis   2.8070451  2.8070451
Skewness  -0.1398503 -0.1398503

$`Statistical inference of z series from iteration`[[2]]
         Unperturbed  Perturbed
Mean      0.01638025 0.01638025
Std Dev   1.33770800 1.33770800
Variance  1.78946270 1.78946270
Kurtosis  2.81149042 2.81149042
Skewness  0.06168373 0.06168373

$`Statistical inference of z series from iteration`[[3]]
         Unperturbed   Perturbed
Mean     -0.03840117 -0.03840117
Std Dev   1.34183108  1.34183108
Variance  1.80051065  1.80051065
Kurtosis  2.54877066  2.54877066
Skewness -0.05682836 -0.05682836

$`Statistical inference of z series from iteration`[[4]]
         Unperturbed   Perturbed
Mean     -0.26398874 -0.26398874
Std Dev   1.34811285  1.34811285
Variance  1.81740825  1.81740825
Kurtosis  2.70324307  2.70324307
Skewness  0.01467297  0.01467297

$`Statistical inference of z series from iteration`[[5]]
         Unperturbed   Perturbed
Mean     -0.05310944 -0.05310944
Std Dev   1.35397535  1.35397535
Variance  1.83324926  1.83324926
Kurtosis  2.56821492  2.56821492
Skewness -0.08533617 -0.08533617

$`Statistical inference of z series from iteration`[[6]]
         Unperturbed   Perturbed
Mean     -0.15963401 -0.15963401
Std Dev   1.29015379  1.29015379
Variance  1.66449680  1.66449680
Kurtosis  2.80113619  2.80113619
Skewness -0.05926602 -0.05926602

$`Statistical inference of z series from iteration`[[7]]
         Unperturbed   Perturbed
Mean     -0.10552435 -0.10552435
Std Dev   1.21692676  1.21692676
Variance  1.48091075  1.48091075
Kurtosis  3.31274071  3.31274071
Skewness  0.03480875  0.03480875

$`Statistical inference of z series from iteration`[[8]]
         Unperturbed   Perturbed
Mean     -0.06228621 -0.06228621
Std Dev   1.35975486  1.35975486
Variance  1.84893329  1.84893329
Kurtosis  2.75713484  2.75713484
Skewness -0.13113248 -0.13113248

$`Statistical inference of z series from iteration`[[9]]
         Unperturbed   Perturbed
Mean     0.141133580 0.141133581
Std Dev  1.399466931 1.399466931
Variance 1.958507692 1.958507692
Kurtosis 2.451595198 2.451595198
Skewness 0.004259234 0.004259234

$`Statistical inference of z series from iteration`[[10]]
         Unperturbed   Perturbed
Mean      0.02717184  0.02717184
Std Dev   1.36086399  1.36086399
Variance  1.85195079  1.85195079
Kurtosis  2.71124354  2.71124354
Skewness -0.15811334 -0.15811334


$`errors of z series from iteration`
$`errors of z series from iteration`[[1]]
       z_errors
1 -9.258905e-10
2 -9.258904e-10
3 -9.258904e-10
4 -9.258904e-10
5 -9.258904e-10
6 -9.258905e-10

$`errors of z series from iteration`[[2]]
       z_errors
1 -7.459174e-10
2 -7.459173e-10
3 -7.459174e-10
4 -7.459174e-10
5 -7.459174e-10
6 -7.459173e-10

$`errors of z series from iteration`[[3]]
       z_errors
1 -9.792326e-10
2 -9.792327e-10
3 -9.792326e-10
4 -9.792326e-10
5 -9.792326e-10
6 -9.792327e-10

$`errors of z series from iteration`[[4]]
       z_errors
1 -1.350503e-11
2 -1.350520e-11
3 -1.350503e-11
4 -1.350520e-11
5 -1.350503e-11
6 -1.350520e-11

$`errors of z series from iteration`[[5]]
       z_errors
1 -5.653198e-10
2 -5.653198e-10
3 -5.653198e-10
4 -5.653198e-10
5 -5.653198e-10
6 -5.653198e-10

$`errors of z series from iteration`[[6]]
       z_errors
1 -6.353736e-10
2 -6.353735e-10
3 -6.353736e-10
4 -6.353736e-10
5 -6.353736e-10
6 -6.353735e-10

$`errors of z series from iteration`[[7]]
       z_errors
1 -9.764009e-10
2 -9.764009e-10
3 -9.764009e-10
4 -9.764008e-10
5 -9.764008e-10
6 -9.764010e-10

$`errors of z series from iteration`[[8]]
       z_errors
1 -4.043135e-10
2 -4.043135e-10
3 -4.043135e-10
4 -4.043135e-10
5 -4.043135e-10
6 -4.043135e-10

$`errors of z series from iteration`[[9]]
       z_errors
1 -3.782752e-10
2 -3.782752e-10
3 -3.782752e-10
4 -3.782752e-10
5 -3.782752e-10
6 -3.782752e-10

$`errors of z series from iteration`[[10]]
       z_errors
1 -1.455423e-10
2 -1.455422e-10
3 -1.455423e-10
4 -1.455424e-10
5 -1.455423e-10
6 -1.455422e-10


$`Statistical inference of errors from iteration`
$`Statistical inference of errors from iteration`[[1]]
                output
Mean     -9.258904e-10
Std Dev   4.819111e-17
Variance  2.322383e-33
Kurtosis  1.340104e+00
Skewness -3.826014e-01

$`Statistical inference of errors from iteration`[[2]]
                output
Mean     -7.459174e-10
Std Dev   7.477871e-17
Variance  5.591855e-33
Kurtosis  2.316277e+00
Skewness -8.406419e-01

$`Statistical inference of errors from iteration`[[3]]
                output
Mean     -9.792326e-10
Std Dev   5.477577e-17
Variance  3.000385e-33
Kurtosis  1.081989e+00
Skewness -2.799286e-01

$`Statistical inference of errors from iteration`[[4]]
                output
Mean     -1.350507e-11
Std Dev   8.381780e-17
Variance  7.025423e-33
Kurtosis  1.848953e+00
Skewness -7.418875e-01

$`Statistical inference of errors from iteration`[[5]]
                output
Mean     -5.653198e-10
Std Dev   2.126461e-17
Variance  4.521835e-34
Kurtosis  1.843398e+00
Skewness -4.707716e-01

$`Statistical inference of errors from iteration`[[6]]
                output
Mean     -6.353736e-10
Std Dev   4.689149e-17
Variance  2.198812e-33
Kurtosis  1.358127e+00
Skewness  3.768862e-01

$`Statistical inference of errors from iteration`[[7]]
                output
Mean     -9.764008e-10
Std Dev   6.416134e-17
Variance  4.116677e-33
Kurtosis  2.942341e+00
Skewness  1.401559e-01

$`Statistical inference of errors from iteration`[[8]]
                output
Mean     -4.043135e-10
Std Dev   2.503721e-17
Variance  6.268619e-34
Kurtosis  1.136152e+00
Skewness  6.245280e-03

$`Statistical inference of errors from iteration`[[9]]
                output
Mean     -3.782752e-10
Std Dev   1.027099e-17
Variance  1.054933e-34
Kurtosis  2.561983e+00
Skewness -1.084713e+00

$`Statistical inference of errors from iteration`[[10]]
                output
Mean     -1.455423e-10
Std Dev   3.479104e-17
Variance  1.210416e-33
Kurtosis  2.236902e+00
Skewness -1.635231e-01
```
So let us explore the output, the first list is the ``z series from iteration``, so for first iteration we have the following z series:
```{coffee}
  Unperturbed   Perturbed
1 -3.41517722 -3.41517722
2  0.98871744  0.98871744
3 -0.03742858 -0.03742858
4 -0.68730508 -0.68730508
5 -0.32776829 -0.32776829
6 -2.43527207 -2.43527207
```
Now, we only have 6 data points here, but actually there are 100 data points, since the defualt number of iterations `r = 100`. We only have 6 here since the default `output` parameter of the function is set to `"partial"`.

So there seems to be something wrong with the z series because even in the first iteration we have the same z series for unperturbed and perturbed, even though we have set the error for perturbed, which is the `epsilon * yt`. And that is where the problem actually, knowing that `epsilon = 10 ^ (-9)`, then if we have `yt =  0.6097608`. Then the perturbed model would be `z + epsilon * yt`, wherein `epsilon * yt = 10^(-9) *  0.6097608 =  6.097608e-10`. So the difference between unperturbed and perturbed would be  `6.097608e-10`, that is so small that we forgot to notice it. That is why in the above output for unperturbed and perturbed model on the first iteration and other iterations the z series are seems to be the same, but there is acutally a very small difference of `10 ^ (-9) * yt`.

Next we took the statistical inference for the z series, and that is the next list of the output, the ``Statistical inference of z series from iteration``. So the statistical inference of z series from the first iteration would be,
```{coffee}
         Unperturbed  Perturbed
Mean      -0.1763497 -0.1763497
Std Dev    1.3492610  1.3492610
Variance   1.8205051  1.8205051
Kurtosis   2.8070451  2.8070451
Skewness  -0.1398503 -0.1398503
```
So what to expect? Since both z series for unperturbed and perturbed have very small difference, then it should be expected that their statistcal inference should be close as well.

The ``errors of z series from iteration`` is just the difference of the between the unperturbed and perturbed z series. Now, since the error of the perturbed model is `epsilon * yt`, where `yt` is generated from the uniform distribution (0, 1). Then for n initial points, we have n errors for perturbed. Of course, for first initial points, we have a corressponding error for perturbed, this error would then be used for `r` iterations, so if `r = 100000`, then there would be 100000 iterations and the first error of the perturbed model, for the first initial points is used for the entire itration before proceeding to the next error on the next intial points. That is why, for the entire iteration we have the same set of errors. From the example, the errors of z series from the first iteration would be,
```{coffee}
1 -9.258905e-10
2 -9.258904e-10
3 -9.258904e-10
4 -9.258904e-10
5 -9.258904e-10
6 -9.258905e-10
```

Also, this error is so obvious without numerically computing it. Since from the theory, the unperturbed model is `z`, and the perturbed model is `z + epsilon * yt`. Hence, the error is just `epsilon * yt`, since that is the difference between the unperturbed and perturbed model.
