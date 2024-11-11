# stat computing123

# 3.2 laplace distribution
# acceptance-rejection method
# f(x) = exp(-|x|)/2
# f(x) < 1/2
# g(x) = 1/12
# c = 6
# determine the upper and lower for x
-optim(0,function (x) -exp(-abs(x))/2)$value
uniroot(function (x) exp(-abs(x))/2 - 0.001,c(-10,0),tol = 1e-9)$root



n <- 1000;k <- j <- 0
y <- numeric(n)

while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1,-6,6)
  if (exp(-abs(x)) > u){
    k <- k+1
    y[k]  <- x
 } 
}
plot(density(y),col='red')
lines(density(VGAM::rlaplace(n)))
# the inverse transform method
# F^{-1}(p)= sgn(p-0.5)ln(1-2|p-0.5|)
set.seed(123)
p = runif(n)
x = -sign(p-0.5)*log(1-2*abs(p-0.5))
plot(density(x),col='red')
lines(density(VGAM::rlaplace(n)))

# exercise 3.3
# F(x) = 1 - (b/x)^a
# x = exp(log(b) - log(1-u)/a)
a <- b <- 2
u <- runif(n)
x = exp(log(b) - log(1-u)/a)
plot(density(x),col='red')
lines(density(VGAM::rpareto(n,2,2)))


# exercise 3.4:the rayleigh density
# F(x)=1-exp(-x^2/(2sigma^2))

library(VGAM)
n = 1000

scale = 1
ans <- scale * sqrt(-2 * log(runif(n)))

plot(density(ans),col='red')
lines(density(rrayleigh(n)))

# 3.5 
set.seed(123)
library(dplyr)
n <- 10^3
u = runif(n)
x = case_when(u <= 0.1 ~ 0,
              u > 0.1 & u <= 0.3 ~ 1,
              u > 0.3 & u <= 0.5 ~ 2,
              u > 0.5 & u <= 0.7 ~ 3,
              u > 0.7 & u <= 1 ~ 4)

y <- sample(0:4, size = n, replace = T, 
       prob = c(0.1,0.2,0.2,0.2,0.3))

rbind(table(x)/n,table(y)/n)

# 3.6 prove

# 3.7
# f(x) < 1.8
# c = 2


myfunc <- function(n) {
k <- j <- 0
y <- numeric(n)

while (k < n) {
  u <- runif(1)
  j <- j + 1
  x <- runif(1)
  if (6*x^2*(1-x) > u){
    k <- k+1
    y[k]  <- x
  } 
}
return(y)
}
plot(density(myfunc(n)),col='red')
lines(density(rbeta(n,3,2)))

# 3.8
rm(list=ls)
myfunc <- function(n,mu,sigma) exp(rnorm(n,mu,sigma))

n <- 10^3;mu <- 0;sigma <- 1

x <- myfunc(n,mu,sigma)
# please note :The log normal distribution has density

# f(x) = 1/(√(2 π) σ x) e^-((log x - μ)^2 / (2 σ^2))

# where μ and σ are the mean and standard deviation of the logarithm. 
# The mean is E(X) = exp(μ + 1/2 σ^2),
# and the variance Var(X) = exp(2*μ + σ^2)*(exp(σ^2) - 1)

plot(x,dlnorm(x,meanlog = mu + sigma^2/2,
              sdlog= 2*mu + sigma^2 - log(exp(sigma^2) - 1)),
     col='black')

lines(density(x),col='red')

# 3.9
rm(list=ls())
myfunc <- function(n){
  U1 <- runif(n,-1,1)
  U2 <- runif(n,-1,1)
  U3 <- runif(n,-1,1)
  return(ifelse(abs(U3) >= abs(U1) & abs(U3) >= abs(U2),U2,U3))
  
}
x <- myfunc(10^3)

plot(x,(1-x^2)*3/4,col='red')
lines(density(x))

# 3.11
rm(list=ls())
p1 <- 0.75;n=1000
y <- p1*rnorm(n) + (1-p1)*rnorm(n,mean = 3,sd = 1)
plot(density(y))

for ( i in seq(0.1,0.9,0.1)){
  lines(density(i*rnorm(n) + (1-i)*rnorm(n,mean = 3,sd = 1)),col= 10*i)
}

# 3.12
n = 10^3;beta=2;r=4
y = rexp(n,rgamma(1,shape = 4,rate = 2))

# 3.13
plot(density(y),col ='red')
lines(density(exp(log(beta)-log(runif(n))/r)-beta))

# 3.14
rm(list=ls())
rmvn.Choleski <-
  function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    Q <- chol(Sigma) # Choleski factorization of Sigma
    Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    X
  }

X <- rmvn.Choleski(200, mu = 0:2, 
                   Sigma = matrix(c(10,-5,5,-5,10,-5,5,-5,10)/10,
                                  nrow=3,byrow = T))
pairs(X)

cor(X)
apply(X,2,mean)


# 3.17
rmvn.Choleski <-
  function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    Q <- chol(Sigma) # Choleski factorization of Sigma
    Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    X
  }
X <- rmvn.Choleski(200, mu, Sigma)
pairs(X)

# 3.18
d = 5
(m2 <- matrix(1:d^2, 5, 5))
# lower.tri(m2)
m2[lower.tri(m2)] <- 0
m2

for ( i in 1:d) {
  for (j in i:d) {
    m2[i,j] = ifelse(i == j,sqrt(rchisq(1,d-i+1)),rnorm(1))
    
  }
}
m2


# Chapter 6
# 6.1
n <- 10^4
x <- runif(n,min = 0,max = pi/3)
theta_hat <- mean(sin(x))*pi/3
print(theta_hat)
print(cos(0)-cos(pi/3))

# 6.2
x <- seq(.1,2.5,length = 10)
m <- 10^4
u <- runif(m)
cdf <- numeric(length(x))
for (i in 1:length(x)){
  cdf[i] <-  mean(x[i]*exp(-(u*x[i])^2/2))/sqrt(2*pi)
  
}
phi <- pnorm(x)-pnorm(0)
print(round(rbind(x,cdf,phi),3))

x <- 2;m <- 10^4;z <- rnorm(m)
g <- (z < x); v <- mean((g - mean(z))^2)/m
cdf <- mean(g)
v;c(cdf - 1.96*sqrt(v),cdf + 1.96*sqrt(v))

# 6.3
set.seed(123)
rm(list = ls())
n <- 10^4;x <- runif(n,min = 0, max = 0.5)

theta.hat <- mean(exp(-x))/2

print(theta.hat);print(1- exp(-0.5))

b  <-  1/2;a  <-  0
var_g <- mean((exp(-x) - mean(exp(-x)))^2)
v <- var_g*(b - a)^2/n
v
# exponential
set.seed(123)
rm(list = ls())
n <- 10^4;x <- rexp(n)

theta.hat <- mean(exp(-x))/2

print(theta.hat);print(1- exp(-0.5))

b  <-  1/2;a  <-  0
var_g <- mean((exp(-x) - mean(exp(-x)))^2)
v <- var_g*(b - a)^2/n
v

# 6.4
rm(list = ls())
set.seed(123)

x <- c(1:9)/10
m <- 10^4
u <- runif(m)
cdf <- numeric(length(x))
for (i in 1:length(x)){
  g <- u^2*x[i]^3*(1 - u*x[i])^2
  cdf[i] <- gamma(6)*mean(g)/gamma(3)^2
}
phi <- pbeta(x,3,3)
print(round(rbind(x,cdf,phi),3))

# 6.5 

rm(list = ls())
set.seed(123)

x <- 2;m <- 10^4;z <- rnorm(m);g <- (z < x)
v <- mean((g-mean(g))^2)/m;v

g1 <- mean(z < x)
v1 <- g1*(1-g1)/m;
v;v1

# 6.6
rm(list = ls())
set.seed(123)
m <- 10^4
u <- runif(m)
cov(exp(u),exp(1-u))
var(exp(u));var(exp(1-u));var(exp(u) + exp(1-u))

# simple MC vs antithetic variates
x <- numeric()
MC.Phi <- function(R = 10^4, antithetic = TRUE) {
    u <- runif(R/2)
    if (!antithetic) v <- runif(R/2) else
      v <- 1 - u
    u <- c(u, v)
    g <- exp(u)
    x <- mean(g) 
  return(x)
}
m <- 10^3
MC1 <- MC2 <- numeric(m)

for (i in 1:m){
  MC1[i] <- MC.Phi(R = 10^3,antithetic = F)
  MC2[i] <- MC.Phi(R = 10^3,antithetic = T)
}

var(MC1)
var(MC2)
(1-var(MC2)/var(MC1))*100  #reduction in variance

# 6.9
rm(list = ls())
Rayleigh <- function(x, m = 10000, antithetic=TRUE){
  u <- runif(m/2)
  if (!antithetic)
    v <- runif(m/2)
  else
    v <- 1 - u
  u <- c(u,v)
  #calculate inv-cdf and varianace
  invcdf <- x * sqrt(-2 * log(u))
  invcdf
}

set.seed(123)
X1 <- Rayleigh(1.95, antithetic=FALSE)
X2 <- Rayleigh(1.95, antithetic=TRUE)
var(X1);var(X2)
# plot(density(X2))
# lines(density(VGAM::rrayleigh(10^4,scale = 1.95)))
set.seed(321)
X1prime <- Rayleigh(1.95, antithetic=FALSE)
X2prime <- Rayleigh(1.95, antithetic=T)
var(X1prime);var(X2prime)

(var(X1) - var(X2))/var(X1)
(var(X1prime) - var(X2prime))/var(X1prime)

# 6.10
# 7.1 Estimate the MSE of the level k trimmed means for random samples
# of size 20 generated from a standard Cauchy distribution. 
# (The target parameter θ is the center or median; the expected value does not exist.)
# Summarize the estimates of MSE in a table for k = 1,2,...,9.

n <- 20
m <- 1000
tmean <- numeric(m)
for (i in 1:m) {
  x <- sort(rcauchy(n))
  tmean[i] <- median(x[2:(n-1)]) # the estimate of theta is median, not mean
}
mse <- mean(tmean^2)
sqrt(sum((tmean- mean(tmean))^2)) / m

rm(list = ls())
n <- 20;K <- n/2- 1;m <- 1000
mse <- matrix(0, n/2, 2)

trimmed.mse <- function(n, m, k) {
  tmedian <- numeric(m)
  for (i in 1:m) {
    x <- sort(rcauchy(n))
    tmedian[i] <- median(x[(k+1):(n-k)]) 
  }
  mse.est <- mean(tmedian^2)
  se.mse <- sqrt(mean((tmedian-mean(tmedian))^2)) / sqrt(m)
  return(c(mse.est, se.mse))
}
for (k in 0:K) {
  mse[k+1, 1:2] <- trimmed.mse(n=n, m=m, k=k)
  }
mse


# 7.2 Plot the empirical power curve for the t-test in Example 7.9, changing
# the alternative hypothesis to H1 : µ /= 500, and keeping the significance
# level α = 0.05.
rm(list = ls())
n <- 20;m <- 1000
mu0 <- 500;sigma <- 100

mu <- c(seq(450, 650, 10)) #alternatives
M <- length(mu)
power <- numeric(M)
for (i in 1:M) {
  mu1 <- mu[i]
  pvalues <- replicate(m, expr = {
    #simulate under alternative mu1
    x <- rnorm(n, mean = mu1, sd = sigma)
    ttest <- t.test(x,
                    alternative = "two.sided", mu = mu0)
    ttest$p.value } )
  power[i] <- mean(pvalues <= .05)
}
se <- sqrt(power * (1-power) / m)

library(ggplot2)
df <- data.frame(mean=mu, power=power,
                 upper=power+2*se, lower=power-2*se)
ggplot(df, aes(x=mean, y=power)) +
  geom_line() +
  geom_vline(xintercept=500, lty=2) +
  geom_hline(yintercept=c(0,.05), lty=1:2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2, lwd=1.5)


# 7.3 Plot the power curves for the t-test in Example 7.9 for sample sizes 10,
# 20, 30, 40, and 50, but omit the standard error bars. Plot the curves
# on the same graph, each in a different color or different line type, and
# include a legend. Comment on the relation between power and sample
# size.
rm(list = ls())
m <- 1000
mu0 <- 500;sigma <- 100

mu <- c(seq(450, 650, 10)) #alternatives
M <- length(mu)

power <- numeric(M);se <- numeric(M);
myfunc <- function(n) {
for (i in 1:M) {
  mu1 <- mu[i]
  pvalues <- replicate(m, expr = {
    #simulate under alternative mu1
    x <- rnorm(n, mean = mu1, sd = sigma)
    ttest <- t.test(x,
                    alternative = "two.sided", mu = mu0)
    ttest$p.value } )
  power[i] <- mean(pvalues <= .05)
    }
return(power)
}

plot(mu,myfunc(10),type = 'l',col = 1)
lines(mu,myfunc(20),type = 'l' ,col = 2)
lines(mu,myfunc(30),type = 'l',col = 3)
lines(mu,myfunc(40),type = 'l',col = 4)
lines(mu,myfunc(50),type = 'l',col = 5)

        
# 7.4 Suppose that X1,...,Xn are a random sample from a lognormal distri
# bution. Construct a 95% confidence interval for the parameter µ. Use a
# Monte Carlo method to obtain an empirical estimate of the confidence
# level when data is generated from standard lognormal.

rm(list = ls())
n <- 20; m <- 1000
alpha <- .05

UCL = numeric(m)
LCL = numeric(m)
for(i in 1:m){
  x <- log(rlnorm(n, meanlog = 0, sdlog = 1))
  cv <- qnorm(1-alpha/2)
  LCL[i] <- mean(x) - cv/sqrt(n)
  UCL[i] <- mean(x) + cv/sqrt(n)
}
sum(LCL < 0 & UCL > 0)
mean(LCL < 0 & UCL > 0)

        
# 7.5  Suppose that X1, . . . , Xn are a random sample from a lognormal distribution. 
# Construct a 95% confidence interval for the parameter µ.
# Use a Monte Carlo method to obtain an empirical estimate of the confidence level when data is generated from standard lognormal.
rm(list = ls());library(purrr)

wei_runs <- replicate(1000, {
  x <- rbinom(1000, size = 1, prob = .5)
  r <- rle(x)
  mytab <- table(r$lengths)
  len <- length(mytab)
  runs <- max(r$lengths)
  w <- as.numeric(mytab[len])
  return(c(runs,w))
}
)

# the probability that the observed maximum run length in [9, 11]
mytab <- table(wei_runs[1,])
mydf <- data.frame(x =  as.numeric(dimnames(mytab)[[1]]),
                      y = as.numeric(mytab)) 
sum(mydf$y[dplyr::between(mydf$x,9,11)]/sum(mydf$y))

(theta <- sum(wei_runs[1,]*wei_runs[2,])/sum(wei_runs[2,]))

sd <- sqrt(mean((wei_runs[1,]-theta)^2*wei_runs[2,]))

c(theta - 1.96*sd, theta + 1.96*sd)」



# 7.6:Suppose a 95% symmetric t-interval is applied to estimate a mean,
# but the sample data are non-normal. Then the probability that the
# confidence interval covers the mean is not necessarily equal to 0.95. Use
# a Monte Carlo experiment to estimate the coverage probability of the
# t-interval for random samples of χ2(2) data with sample size n = 20.
# Compare your t-interval results with the simulation results in Example
# 7.4. (The t-interval should be more robust to departures from normality
#       than the interval for variance.)

alpha = 0.05
n = 20
m = 1000

UCL = numeric(m)
LCL = numeric(m)

for(i in 1:m)
{
  x = rchisq(n, 2) # compare with x = rnorm(n) + 2
  LCL[i] = mean(x) - qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
  UCL[i] = mean(x) + qt(alpha / 2, df=n-1, lower.tail = FALSE)*sd(x)/sqrt(n)
}

mean(LCL < 2 & UCL > 2)



# 7.7:Estimate the 0.025, 0.05, 0.95, and 0.975 quantiles of the skewness
# √b1 under normality by a Monte Carlo experiment. 
# Compute the standard error of the estimates from (2.14) 
# using the normal approximation for the density (with exact variance formula). 
# Compare the estimated quantiles with the quantiles of the large sample approximation
# √b1 ≈ N(0, 6/n).

rm(list = ls())
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

n <- 20;m <- 1000
y <- numeric(m)
for (i in 1:m) {
  x <- rnorm(n)
  y[i] <- sk(x)
}

quantile(y,probs = c(.025,.05,.95,.975))
qnorm(c(.025,.05,.95,.975),mean = 0, sd = sqrt(6*(n-2)/((n+1)*(n+3))))
qnorm(c(.025,.05,.95,.975),mean = 0, sd = sqrt(6/n))



# 7.8：Estimate the power of the skewness test of normality against 
# symmetric Beta(α, α) distributions and comment on the results. 
# Are the results different for heavy-tailed symmetric alternatives such as t(ν)?

alpha <- .1
n <- 30
m <- 2500

#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

mean(replicate(m, expr = {
  #simulate under alternative mu1
  x <- rnorm(n, mean = 0, sd = 1)
  as.integer(abs(sk(x)) >= cv)} ))

mean(replicate(m, expr = {
  #simulate under alternative mu1
  x <- rbeta(n, 0.5, 0.5)
  as.integer(abs(sk(x)) >= cv)} ))

mean(replicate(m, expr = {
  #simulate under alternative mu1
  x <- rt(n, 10)
  as.integer(abs(sk(x)) >= cv)} ))


# 7.9
# Refer to Example 7.16. Repeat the simulation, but also compute the
# F test of equal variance, at significance level αˆ = 0.055. Compare the
# power of the Count Five test and F test for small, medium, and large sample sizes. 
# (Recall that the F test is not applicable for non-normal distributions.)

# generate samples under H1 to estimate power
rm(list = ls())
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
sigma1 <- 1
sigma2 <- 1.5
m <- 1000;n <- c(1,10,100)*20
myfunc <- function(i){
print(c(mean(replicate(m, expr={
  x <- rnorm(n[i], 0, sigma1)
  y <- rnorm(n[i], 0, sigma2)
  count5test(x, y)
})),
mean(replicate(m, expr={
  x <- rnorm(n[i], 0, sigma1)
  y <- rnorm(n[i], 0, sigma2)
  var.test(x, y,conf.level = 1- 0.055)$p.value
}))))
}
options(digits = 3)
myfunc(1)
myfunc(2)
myfunc(3)




