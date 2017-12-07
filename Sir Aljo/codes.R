# simulate poisson distribution with mean \lambda

# The skewness of the poisson distribution is
# \frac{1}{\sqrt{\lambda}}, hence a function of \lambda.
# So that in terms of the skewness we have \lambda = \frac{1}{\skewness^2}.

# So if skewness is 3, then \lambda = \frac{1}{3^2} = \frac{1}{9}.
# Thus to simulate poisson with skewness equal to 3 of size 100, we have
hist(rpois(1000, 1/9), breaks = 100)
