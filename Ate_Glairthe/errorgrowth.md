## Statistical Inference
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
## Error Growth
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
    z_errors_list[[i]] <- z_errors
    z_errors_moments[[i]] <- stat_moments(z_errors)
  }
  
  if (output == "partial") {
    z_vec_list <- lapply(z_vec_list, function(x) head(x, n = 10))
    z_errors_moments <- lapply(z_errors_moments, function(x) head(x, n = 10))
  } else if (output == "full") {
    z_vec_list <- z_vec_list; z_errors_moments <- z_errors_moments
  }
  
  return(list(z_list = z_vec_list, z_moments = z_vec_moments, z_errors = z_errors_list,
         z_errors_moments = z_errors_moments))
}
```
