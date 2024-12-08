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

        
# 7.10
n <- 1000
y <- replicate(1000,expr = {
x <- sort(rlnorm(n,meanlog = 0,sdlog = 1))
mu <- mean(x)
sum((2*(1:n) -n -1)*x)/(n^2*mu)
})
plot(density(y))
mean(y);median(y);quantile(y,probs = (1:10)/10)



y_unif <- replicate(1000,expr = {
  x <- sort(runif(n))
  mu <- mean(x)
  sum((2*(1:n) -n -1)*x)/(n^2*mu)
})

plot(density(y_unif))


y_binom <- replicate(1000,expr = {
  x <- sort(rbinom(n,size = 1,prob = 0.1))
  mu <- mean(x)
  sum((2*(1:n) -n -1)*x)/(n^2*mu)
})

plot(density(y_binom))

# 7.11
rm(list = ls())
n <- 20;alpha <- .05

p_chisq <- replicate(1000,expr = {
  x <- rchisq(n, df = 1)
  t.test(x, alternative = "two.sided", mu = 1)$p.value
})
p_unif <- replicate(1000,expr = {
  x <- runif(n, min = 0,max = 2)
  t.test(x, alternative = "two.sided", mu = 1)$p.value
})
p_exp <- replicate(1000,expr = {
  x <- runif(n, min = 0,max = 2)
  t.test(x, alternative = "two.sided", mu = 1)$p.value
})

mean(p_chisq < alpha);mean(p_unif < alpha);mean(p_exp < alpha)

# 7.12
# two vars are independent
alpha <- .05
cor_pearson <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = diag(c(1,1)) )
  cor.test(x[,1],x[,2],method = "pearson")$p.value
  })

cor_kendall <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = diag(c(1,1)) )
  cor.test(x[,1],x[,2],method = "kendall")$p.value
})

cor_spearman <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = diag(c(1,1)) )
  cor.test(x[,1],x[,2],method = "spearman")$p.value
})

mean(cor_pearson < alpha);
mean(cor_kendall < alpha);
mean(cor_spearman < alpha)
# two vars are dependent
sigma <- matrix(c(1,0.5,0.5,1),nrow = 2)
cor_pearson2 <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = sigma )
  cor.test(x[,1],x[,2],method = "pearson")$p.value
})

cor_kendall2 <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = sigma )
  cor.test(x[,1],x[,2],method = "kendall")$p.value
})

cor_spearman2 <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = sigma )
  cor.test(x[,1],x[,2],method = "spearman")$p.value
})

mean(cor_pearson2 < alpha);
mean(cor_kendall2 < alpha);
mean(cor_spearman2 < alpha)

# 8.1 Compute a jackknife estimate of the bias and the standard error of the
# # correlation statistic in Example 8.2.
rm(list = ls())
# set.seed(123)
myfunc <- function(n,y,z){
  #compute the jackknife replicates, leave-one-out estimates
  theta.jack <- numeric(n)
  for (i in 1:n) theta.jack[i] <- cor(y[-i],z[-i])
  bias <- (n - 1) * (mean(theta.jack) - theta.hat)
  se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2))
  return(c(bias,se))
}
data(patch, package = "bootstrap")
n <- nrow(patch)
y <- patch$y
z <- patch$z
theta.hat <- cor(y,z)

myfunc(n = n, y = y, z = z)

# 
# 8.2 Refer to the law data (bootstrap). Use the jackknife-after-bootstrap
# method to estimate the standard error of the bootstrap estimate of se(R).
data(law,package = 'bootstrap')
names(law)
myfunc(n = n, y = law$LSAT, z = law$GPA)

# 8.3 Obtain a bootstrap t confidence interval estimate for the correlation
# statistic in Example 8.2 (law data in bootstrap).
rm(list = ls())
boot.t.ci <-
  function(x, B = 500, R = 100, level = .95, statistic){
    #compute the bootstrap t CI
    x <- as.matrix(x); n <- nrow(x)
    stat <- numeric(B); se <- numeric(B)
    
    boot.se <- function(x, R, f) {
      #local function to compute the bootstrap
      #estimate of standard error for statistic f(x)
      x <- as.matrix(x); m <- nrow(x)
      th <- replicate(R, expr = {
        i <- sample(1:m, size = m, replace = TRUE)
        f(x[i, ])
      })
      return(sd(th))
    }
    for (b in 1:B) {
      j <- sample(1:n, size = n, replace = TRUE)
      y <- x[j, ]
      stat[b] <- statistic(y)
      se[b] <- boot.se(y, R = R, f = statistic)
    }
    stat0 <- statistic(x)
    t.stats <- (stat - stat0) / se
    se0 <- sd(stat)
    alpha <- 1 - level
    Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
    names(Qt) <- rev(names(Qt))
    CI <- rev(stat0 - Qt * se0)
  }

data(law,package = 'bootstrap')
dat <- cbind(law$LSAT,law$GPA)
stat <- function(dat) {
  cor(dat[, 1],dat[, 2]) }
ci <- boot.t.ci(dat, statistic = stat, B=2000, R=200)
print(ci)


# 8.4 Refer to the air-conditioning data set aircondit provided in the boot package. The 12 observations are the times in hours between failures of
# air-conditioning equipment [68, Example 1.1]:
#   3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487.
# Assume that the times between failures follow an exponential model
# Exp(λ). Obtain the MLE of the hazard rate λ and use bootstrap to
# estimate the bias and standard error of the estimate.
rm(list = ls())
data(aircondit,package = 'boot')
time <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
lambda <- 12/mean(time)

#set up the bootstrap
B <- 200 #number of replicates
n <- length(time) #sample size
R <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  lambda <- time[i] #i is a vector of indices
  R[b] <- 12/mean(lambda)
}
#output
print(se.R <- sd(R))
hist(R, prob = TRUE)



rm(list = ls())
data(aircondit,package = 'boot')
time <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
inv_lambda <- mean(time)/12

# 8.5 Refer to Exercise 8.4. Compute 95% bootstrap confidence intervals for
# the mean time between failures 1/λ by the standard normal, basic, percentile, and BCa methods. 
# Compare the intervals and explain why they may differ.
rm(list = ls())
data(aircondit,package = 'boot')
# type:c("norm","basic", "stud", "perc", "bca")
boot.out <- boot(aircondit,R = 2000,statistic = function(x,i) mean(x[i,1]))
library(boot)
boot.ci(boot.out, conf = 0.95, 
        type = c("norm","basic", "stud", "perc", "bca"))


# 8.6 Efron and Tibshirani discuss the scor (bootstrap) test score data
# on 88 students who took examinations in five subjects [91, Table 7.1],
# [194, Table 1.2.1]. The first two tests (mechanics, vectors) were closed
# book and the last three tests (algebra, analysis, statistics) were open
# book. Each row of the data frame is a set of scores (xi1, . . . , xi5) for
# the ith student. Use a panel display to display the scatter plots for each pair of test scores. 
# Compare the plot with the sample correlation matrix. 
# Obtain bootstrap estimates of the standard errors for each of the following estimates: 
#   ρˆ12 = ˆρ(mec, vec), ρˆ34 = ˆρ(alg, ana), ρˆ35 = ˆρ(alg,sta), ρˆ45 = ˆρ(ana, sta).

library(bootstrap)
data(scor)

car::scatterplotMatrix(scor)
cor(scor)

#set up the bootstrap
B <- 200 #number of replicates
n <- nrow(scor) #sample size
R12 <- R34 <- R35 <- R45 <- numeric(B) #storage for replicates

#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  mec <- scor$mec[i]
  vec <- scor$vec[i]
  alg <- scor$alg[i]
  ana <- scor$ana[i]
  sta  <- scor$sta[i]
  
  R12[b] <- cor(mec,vec)
  R34[b] <- cor(alg,ana)
  R35[b] <- cor(alg,sta)
  R45[b] <- cor(ana,sta)
}
sd(R12);sd(R34);sd(R35);sd(R45)


# 8.7 Refer to Exercise 8.6. Efron and Tibshirani discuss the following example [91, Chapter 7]. The five-dimensional scores data have a 5 × 5
# covariance matrix Σ, with positive eigenvalues λ1 > · · · > λ5. In principal components analysis,
# measures the proportion of variance explained by the first principal
# component. Let ˆλ1 > · · · > ˆλ5 be the eigenvalues of ˆΣ, where ˆΣ is the
# MLE of Σ. Compute the sample estimate of θ. 
# Use bootstrap to estimate the bias and standard error of θˆ.
#set up the bootstrap
B <- 200 #number of replicates
n <- nrow(scor) #sample size
eg <- eigen(cov(scor))$values
theta_hat <- max(eg)/sum(eg)
theta_b <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  eg <- eigen(cov(scor[i,]))$values
  theta_b[b] <- max(eg)/sum(eg)
}

(bias <- mean(theta_b) - theta_hat)
(se <- sd(theta_b))



# 8.8 Refer to Exercise 8.7. Obtain the jackknife estimates of bias and standard error of θ

B <- 200 #number of replicates
n <- nrow(scor) #sample size
eg <- eigen(cov(scor))$values
theta_hat <- max(eg)/sum(eg)
theta_jack <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
for (b in 1:B) {
  eg <- eigen(cov(scor[-b,]))$values
  theta_jack[b] <- max(eg)/sum(eg)
}

(bias <- (n-1)*(mean(theta_jack) - theta_hat))

(se <- sqrt((n-1)*mean((theta_jack - mean(theta_jack))^2)))


# 8.9 Refer to Exercise 8.7. Compute 95% percentile and BCa confidence
# intervals for θˆ.
theta.boot <- function(dat,ind){
  eg <- eigen(cov(dat[ind,]))$values
  theta_hat <- max(eg)/sum(eg)
}
dat <- scor
boot.out <- boot(dat,statistic = theta.boot,R = 2000)
library(boot)
boot.ci(boot.out, type = c("perc", "bca"))

# 8.10 In Example 8.17, leave-one-out (n-fold) cross validation was used to
# select the best fitting model. Repeat the analysis replacing the Log-Log model with a cubic polynomial model. Which of the four models
# is selected by the cross validation procedure? Which model is selected
# according to maximum adjusted R2 ?

data(ironslag,package = 'DAAG')
names(ironslag)
chemical <- ironslag$chemical
magnetic <- ironslag$magnetic

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  e1[k] <- summary(J1)$adj.r.squared
  J2 <- lm(y ~ x + I(x^2))
  e2[k] <- summary(J2)$adj.r.squared
  J3 <- lm(log(y) ~ x)
  e3[k] <- summary(J3)$adj.r.squared
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  e4[k] <- summary(J4)$adj.r.squared
}
c(mean(e1), mean(e2), mean(e3), mean(e4))

# 8.11 In Example 8.17, leave-one-out (n-fold) cross validation was used to select the best fitting model. 
# Use leave-two-out cross validation to compare the models.


data(ironslag,package = 'DAAG')
names(ironslag)
chemical <- ironslag$chemical
magnetic <- ironslag$magnetic
# 集合运算的一般规则如下：
# union(x,y)    #求并集
# intersect(x,y)    #求交集
# setdiff(x,y)    #求属于x而不属于y的所有元素
# setequal(x,y)    #判断x与y是否相等
# a %in% y    #判断a是否为y中的元素
# choose(n, k)    #n个里面取k个的组合数
# combn(x,n)    #x中的元素每次取n个的所有组合
# combn(x,n,f)     #将这些组合用于指定函数f

n <- 200
adjr1 <- combn(n, 2, function(k) summary(lm(magnetic[-k] ~  chemical[-k]))$adj.r.squared) 
adjr2 <- combn(n, 2, function(k) summary(lm(magnetic[-k] ~  chemical[-k] + I(chemical[-k]^2)))$adj.r.squared) 
adjr3 <- combn(n, 2, function(k) summary(lm(log(magnetic[-k]) ~  chemical[-k]))$adj.r.squared) 
adjr4 <- combn(n, 2, function(k) summary(lm(log(magnetic[-k]) ~  log(chemical[-k])))$adj.r.squared) 

c(mean(adjr1),mean(adjr2),mean(adjr3),mean(adjr4))

# 8.A Conduct a Monte Carlo study to estimate the coverage probabilities of
# the standard normal bootstrap confidence interval, the basic bootstrap
# confidence interval, and the percentile confidence interval. Sample from
# a normal population and check the empirical coverage rates for the
# sample mean. Find the proportion of times that the confidence intervals
# miss on the left, and the porportion of times that the confidence intervals
# miss on the right.
# 8.B Repeat Project A for the sample skewness statistic. 
# Compare the coverage rates for normal populations (skewness 0) and χ2 (5) distributions (positive skewness).

library(boot)
data(law, package = "bootstrap")
boot.obj <- boot(law, R = 2000,
                 statistic = function(x, i){cor(x[i,1], x[i,2])})
print(boot.ci(boot.obj, type=c("basic","norm","perc")))

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

        
# 7.10
n <- 1000
y <- replicate(1000,expr = {
x <- sort(rlnorm(n,meanlog = 0,sdlog = 1))
mu <- mean(x)
sum((2*(1:n) -n -1)*x)/(n^2*mu)
})
plot(density(y))
mean(y);median(y);quantile(y,probs = (1:10)/10)



y_unif <- replicate(1000,expr = {
  x <- sort(runif(n))
  mu <- mean(x)
  sum((2*(1:n) -n -1)*x)/(n^2*mu)
})

plot(density(y_unif))


y_binom <- replicate(1000,expr = {
  x <- sort(rbinom(n,size = 1,prob = 0.1))
  mu <- mean(x)
  sum((2*(1:n) -n -1)*x)/(n^2*mu)
})

plot(density(y_binom))

# 7.11
rm(list = ls())
n <- 20;alpha <- .05

p_chisq <- replicate(1000,expr = {
  x <- rchisq(n, df = 1)
  t.test(x, alternative = "two.sided", mu = 1)$p.value
})
p_unif <- replicate(1000,expr = {
  x <- runif(n, min = 0,max = 2)
  t.test(x, alternative = "two.sided", mu = 1)$p.value
})
p_exp <- replicate(1000,expr = {
  x <- runif(n, min = 0,max = 2)
  t.test(x, alternative = "two.sided", mu = 1)$p.value
})

mean(p_chisq < alpha);mean(p_unif < alpha);mean(p_exp < alpha)

# 7.12
# two vars are independent
alpha <- .05
cor_pearson <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = diag(c(1,1)) )
  cor.test(x[,1],x[,2],method = "pearson")$p.value
  })

cor_kendall <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = diag(c(1,1)) )
  cor.test(x[,1],x[,2],method = "kendall")$p.value
})

cor_spearman <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = diag(c(1,1)) )
  cor.test(x[,1],x[,2],method = "spearman")$p.value
})

mean(cor_pearson < alpha);
mean(cor_kendall < alpha);
mean(cor_spearman < alpha)
# two vars are dependent
sigma <- matrix(c(1,0.5,0.5,1),nrow = 2)
cor_pearson2 <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = sigma )
  cor.test(x[,1],x[,2],method = "pearson")$p.value
})

cor_kendall2 <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = sigma )
  cor.test(x[,1],x[,2],method = "kendall")$p.value
})

cor_spearman2 <- replicate(1000,expr = {
  x <- MASS::mvrnorm(n, mu = c(0,0),Sigma = sigma )
  cor.test(x[,1],x[,2],method = "spearman")$p.value
})

mean(cor_pearson2 < alpha);
mean(cor_kendall2 < alpha);
mean(cor_spearman2 < alpha)

# 8.1 Compute a jackknife estimate of the bias and the standard error of the
# # correlation statistic in Example 8.2.
rm(list = ls())
# set.seed(123)
myfunc <- function(n,y,z){
  #compute the jackknife replicates, leave-one-out estimates
  theta.jack <- numeric(n)
  for (i in 1:n) theta.jack[i] <- cor(y[-i],z[-i])
  bias <- (n - 1) * (mean(theta.jack) - theta.hat)
  se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2))
  return(c(bias,se))
}
data(patch, package = "bootstrap")
n <- nrow(patch)
y <- patch$y
z <- patch$z
theta.hat <- cor(y,z)

myfunc(n = n, y = y, z = z)

# 
# 8.2 Refer to the law data (bootstrap). Use the jackknife-after-bootstrap
# method to estimate the standard error of the bootstrap estimate of se(R).
data(law,package = 'bootstrap')
names(law)
myfunc(n = n, y = law$LSAT, z = law$GPA)

# 8.3 Obtain a bootstrap t confidence interval estimate for the correlation
# statistic in Example 8.2 (law data in bootstrap).
rm(list = ls())
boot.t.ci <-
  function(x, B = 500, R = 100, level = .95, statistic){
    #compute the bootstrap t CI
    x <- as.matrix(x); n <- nrow(x)
    stat <- numeric(B); se <- numeric(B)
    
    boot.se <- function(x, R, f) {
      #local function to compute the bootstrap
      #estimate of standard error for statistic f(x)
      x <- as.matrix(x); m <- nrow(x)
      th <- replicate(R, expr = {
        i <- sample(1:m, size = m, replace = TRUE)
        f(x[i, ])
      })
      return(sd(th))
    }
    for (b in 1:B) {
      j <- sample(1:n, size = n, replace = TRUE)
      y <- x[j, ]
      stat[b] <- statistic(y)
      se[b] <- boot.se(y, R = R, f = statistic)
    }
    stat0 <- statistic(x)
    t.stats <- (stat - stat0) / se
    se0 <- sd(stat)
    alpha <- 1 - level
    Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
    names(Qt) <- rev(names(Qt))
    CI <- rev(stat0 - Qt * se0)
  }

data(law,package = 'bootstrap')
dat <- cbind(law$LSAT,law$GPA)
stat <- function(dat) {
  cor(dat[, 1],dat[, 2]) }
ci <- boot.t.ci(dat, statistic = stat, B=2000, R=200)
print(ci)


# 8.4 Refer to the air-conditioning data set aircondit provided in the boot package. The 12 observations are the times in hours between failures of
# air-conditioning equipment [68, Example 1.1]:
#   3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487.
# Assume that the times between failures follow an exponential model
# Exp(λ). Obtain the MLE of the hazard rate λ and use bootstrap to
# estimate the bias and standard error of the estimate.
rm(list = ls())
data(aircondit,package = 'boot')
time <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
lambda <- 12/mean(time)

#set up the bootstrap
B <- 200 #number of replicates
n <- length(time) #sample size
R <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  lambda <- time[i] #i is a vector of indices
  R[b] <- 12/mean(lambda)
}
#output
print(se.R <- sd(R))
hist(R, prob = TRUE)



rm(list = ls())
data(aircondit,package = 'boot')
time <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
inv_lambda <- mean(time)/12

# 8.5 Refer to Exercise 8.4. Compute 95% bootstrap confidence intervals for
# the mean time between failures 1/λ by the standard normal, basic, percentile, and BCa methods. 
# Compare the intervals and explain why they may differ.
rm(list = ls())
data(aircondit,package = 'boot')
# type:c("norm","basic", "stud", "perc", "bca")
boot.out <- boot(aircondit,R = 2000,statistic = function(x,i) mean(x[i,1]))
library(boot)
boot.ci(boot.out, conf = 0.95, 
        type = c("norm","basic", "stud", "perc", "bca"))


# 8.6 Efron and Tibshirani discuss the scor (bootstrap) test score data
# on 88 students who took examinations in five subjects [91, Table 7.1],
# [194, Table 1.2.1]. The first two tests (mechanics, vectors) were closed
# book and the last three tests (algebra, analysis, statistics) were open
# book. Each row of the data frame is a set of scores (xi1, . . . , xi5) for
# the ith student. Use a panel display to display the scatter plots for each pair of test scores. 
# Compare the plot with the sample correlation matrix. 
# Obtain bootstrap estimates of the standard errors for each of the following estimates: 
#   ρˆ12 = ˆρ(mec, vec), ρˆ34 = ˆρ(alg, ana), ρˆ35 = ˆρ(alg,sta), ρˆ45 = ˆρ(ana, sta).

library(bootstrap)
data(scor)

car::scatterplotMatrix(scor)
cor(scor)

#set up the bootstrap
B <- 200 #number of replicates
n <- nrow(scor) #sample size
R12 <- R34 <- R35 <- R45 <- numeric(B) #storage for replicates

#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  mec <- scor$mec[i]
  vec <- scor$vec[i]
  alg <- scor$alg[i]
  ana <- scor$ana[i]
  sta  <- scor$sta[i]
  
  R12[b] <- cor(mec,vec)
  R34[b] <- cor(alg,ana)
  R35[b] <- cor(alg,sta)
  R45[b] <- cor(ana,sta)
}
sd(R12);sd(R34);sd(R35);sd(R45)


# 8.7 Refer to Exercise 8.6. Efron and Tibshirani discuss the following example [91, Chapter 7]. The five-dimensional scores data have a 5 × 5
# covariance matrix Σ, with positive eigenvalues λ1 > · · · > λ5. In principal components analysis,
# measures the proportion of variance explained by the first principal
# component. Let ˆλ1 > · · · > ˆλ5 be the eigenvalues of ˆΣ, where ˆΣ is the
# MLE of Σ. Compute the sample estimate of θ. 
# Use bootstrap to estimate the bias and standard error of θˆ.
#set up the bootstrap
B <- 200 #number of replicates
n <- nrow(scor) #sample size
eg <- eigen(cov(scor))$values
theta_hat <- max(eg)/sum(eg)
theta_b <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  eg <- eigen(cov(scor[i,]))$values
  theta_b[b] <- max(eg)/sum(eg)
}

(bias <- mean(theta_b) - theta_hat)
(se <- sd(theta_b))



# 8.8 Refer to Exercise 8.7. Obtain the jackknife estimates of bias and standard error of θ

B <- 200 #number of replicates
n <- nrow(scor) #sample size
eg <- eigen(cov(scor))$values
theta_hat <- max(eg)/sum(eg)
theta_jack <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
for (b in 1:B) {
  eg <- eigen(cov(scor[-b,]))$values
  theta_jack[b] <- max(eg)/sum(eg)
}

(bias <- (n-1)*(mean(theta_jack) - theta_hat))

(se <- sqrt((n-1)*mean((theta_jack - mean(theta_jack))^2)))


# 8.9 Refer to Exercise 8.7. Compute 95% percentile and BCa confidence
# intervals for θˆ.
theta.boot <- function(dat,ind){
  eg <- eigen(cov(dat[ind,]))$values
  theta_hat <- max(eg)/sum(eg)
}
dat <- scor
boot.out <- boot(dat,statistic = theta.boot,R = 2000)
library(boot)
boot.ci(boot.out, type = c("perc", "bca"))

# 8.10 In Example 8.17, leave-one-out (n-fold) cross validation was used to
# select the best fitting model. Repeat the analysis replacing the Log-Log model with a cubic polynomial model. Which of the four models
# is selected by the cross validation procedure? Which model is selected
# according to maximum adjusted R2 ?

data(ironslag,package = 'DAAG')
names(ironslag)
chemical <- ironslag$chemical
magnetic <- ironslag$magnetic

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  e1[k] <- summary(J1)$adj.r.squared
  J2 <- lm(y ~ x + I(x^2))
  e2[k] <- summary(J2)$adj.r.squared
  J3 <- lm(log(y) ~ x)
  e3[k] <- summary(J3)$adj.r.squared
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  e4[k] <- summary(J4)$adj.r.squared
}
c(mean(e1), mean(e2), mean(e3), mean(e4))

# 8.11 In Example 8.17, leave-one-out (n-fold) cross validation was used to select the best fitting model. 
# Use leave-two-out cross validation to compare the models.


data(ironslag,package = 'DAAG')
names(ironslag)
chemical <- ironslag$chemical
magnetic <- ironslag$magnetic
# 集合运算的一般规则如下：
# union(x,y)    #求并集
# intersect(x,y)    #求交集
# setdiff(x,y)    #求属于x而不属于y的所有元素
# setequal(x,y)    #判断x与y是否相等
# a %in% y    #判断a是否为y中的元素
# choose(n, k)    #n个里面取k个的组合数
# combn(x,n)    #x中的元素每次取n个的所有组合
# combn(x,n,f)     #将这些组合用于指定函数f

n <- 200
adjr1 <- combn(n, 2, function(k) summary(lm(magnetic[-k] ~  chemical[-k]))$adj.r.squared) 
adjr2 <- combn(n, 2, function(k) summary(lm(magnetic[-k] ~  chemical[-k] + I(chemical[-k]^2)))$adj.r.squared) 
adjr3 <- combn(n, 2, function(k) summary(lm(log(magnetic[-k]) ~  chemical[-k]))$adj.r.squared) 
adjr4 <- combn(n, 2, function(k) summary(lm(log(magnetic[-k]) ~  log(chemical[-k])))$adj.r.squared) 

c(mean(adjr1),mean(adjr2),mean(adjr3),mean(adjr4))
# the other method


# 8.A Conduct a Monte Carlo study to estimate the coverage probabilities of
# the standard normal bootstrap confidence interval, the basic bootstrap
# confidence interval, and the percentile confidence interval. Sample from
# a normal population and check the empirical coverage rates for the
# sample mean. Find the proportion of times that the confidence intervals
# miss on the left, and the porportion of times that the confidence intervals
# miss on the right.

# 8.B Repeat Project A for the sample skewness statistic. 
# Compare the coverage rates for normal populations (skewness 0) and χ2 (5) distributions (positive skewness).

library(boot)
data(law, package = "bootstrap")
boot.obj <- boot(law, R = 2000,
                 statistic = function(x, i){cor(x[i,1], x[i,2])})
print(boot.ci(boot.obj, type=c("basic","norm","perc")))

# 9.1 In jackknife-after-bootstrap, show that the probability that a bootstrap
# sample omits a given observation i is asymptotically equal to e^ −1
# (about 0.368). That is, show that for large n the proportion of bootstrap samples 
# that omit observation i is approximately equal to e^−1, i = 1, . . . , n.
library(boot)
library(bootstrap)
set.seed(1111)
theta.boot <- function(patch, i) {
  # function to compute the patch ratio statistic
  y <- patch[i, "y"]
  z <- patch[i, "z"]
  mean(y) / mean(z)
}
boot.out <- boot(bootstrap::patch,
                 statistic = theta.boot, R=2000)
A <- boot.array(boot.out)
mean(A[, 1] == 0)

# 9.2 
# a. Fit the simple linear model to the DAAG::ironslag data discussed in Example 8.16.
library(DAAG)
data(ironslag)
L <- lm(magnetic ~ chemical,data = ironslag)
x <- ironslag$chemical
summary(L)
# b. Compute the modified residuals two ways and show that they are equal. 
# Hint: Try using the all.equal function to compare the results.
# i. Method 1: Use the definition and the hatvalues function.
# ii. Method 2: Use the rstandard function with sd = 1.
r1 <- L$residuals/sqrt(1 - hatvalues(L)) #Method 1
r2 <- rstandard(L,sd = 1)  #Method 2
all.equal(r,r2)
# c. Compute the hat values directly from the formula for hjj . 
# Check that your vector h is identical to the hatvalues function result.
L.hat <- hatvalues(L)
x <- ironslag$chemical
y <- ironslag$magnetic
ssx <- (length(x)-1)*var(x)
hjj <- 1/length(x) + (x-mean(x))^2/ssx
all.equal(L.hat,hjj)



# 9.3 Refer to the catsM data in the boot package.
data(catsM,package = 'boot')
# a. Display a fitted line plot (use basic R graphics) for the simple linear regression model 
# predicting body weight (Bwt) from heart weight (Hwt).
attach(catsM)
names(catsM)
plot(Hwt,Bwt)
# b. Display the fitted line plot using ggplot2.
library(ggplot2)
ggplot(data = catsM, mapping = aes(x = Hwt, y = Bwt)) +
  geom_point() + geom_smooth(method = "lm")

# c. Display a plot of residuals vs. fits.
fit <- lm(Bwt ~ Hwt)
plot(fit$fitted.values,fit$residuals)
abline(h=0)

# d. Comment on the fit of this model. Are there any outliers? 
# If so,identify these points by observation number.
# Method 1: Use the Interquartile Range
#find Q1, Q3, and interquartile range for values in points column
Q1 <- quantile(catsM$Bwt, .25)
Q3 <- quantile(catsM$Bwt, .75)
IQR <- IQR(catsM$Bwt)

#subset data where points value is outside 1.5*IQR of Q1 and Q3
outliers <- subset(catsM, catsM$Bwt < (Q1 - 1.5*IQR) | catsM$Bwt > (Q3 + 1.5*IQR))

# Method 2: Use Z-Scores
#create new column that calculates z-score of each value in points column
catsM$z <- scale(catsM$Bwt)
#subset data frame where z-score of points value is greater than 3
outliers <- catsM[catsM$z > 3, ]


# e. Based on your analysis above, to analyze the fit using bootstrap, 
# choose a resampling method and explain your reasoning.
library(boot)
m <- 2000
stats <- function(dat, i) {
  x <- dat$Hwt[i]
  y <- dat$Bwt[i]
  Lb <- lm(y ~ x)
  s <- summary(Lb)$sigma
  c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
}
boot.out <- boot(catsM, statistic = stats, R = 2000)
boot.out$t0

# f. Bootstrap the slopes of this model and 
# obtain a bootstrap estimate the standard error of βˆ1.
boot.out

# g. Use jackknife-after-bootstrap to identify influential observations.
library(boot)
stats <- function(dat, i) {
  x <- dat$Hwt[i]
  y <- dat$Bwt[i]
  Lb <- lm(y ~ x)
  slope = Lb$coeff[2]
  return(slope)
}
boot.out <- boot(catsM, stats, R = 2000)
jack.after.boot(boot.out, useJ=TRUE, stinf=FALSE)

# 9.4 Implement the resampling cases method on the MASS::mammals data using the boot function. 
# Compare your results for the bias and se.
names(mammals)
data(mammals,package ='MASS')
x <- mammals$body
y <- mammals$brain
m <- 2000
n <- NROW(x)
L1 <- lm(y ~ x) #estimate the model
b0 <- L1$coeff[1]; b1 <- L1$coeff[2]
## run bootstrap of cases
out <- replicate(m, expr={
  i <- sample(1:n, replace=TRUE, size=n)
  xstar <- x[i]
  ystar <- y[i]
  Lb <- lm(ystar ~ xstar)
  s <- summary(Lb)$sigma
  c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
})
bootCases <- t(out)
meanCases <- colMeans(bootCases)
sdCases <- apply(bootCases, 2, "sd")

biasInt <- mean(bootCases[,1] - b0) #bias for intercept
biasSlope <- mean(bootCases[,2] - b1) #bias for slope

# 9.5 To investigate the error distribution and influence,
# we would like to have bootstrapped estimates of the squared error. 
# Suppose that a data set contains the response y and predictor x. 
# Write a statistic function for use with the boot function that will return the MSE
# for the fitted simple linear regression model y = β0 + β1x + ε. 
# Test this function by generating random bivariate normal data and running an ordinary
# bootstrap of the MSE for the regression model.

## run bootstrap of cases
rm(list = ls());set.seed(1)
m <- 2000
bi.norm <- matrix(rnorm(4000), 2000, 2)
x <- bi.norm[,1];y <- bi.norm[,2]
n <- NROW(x)

out <- replicate(m, expr={
  i <- sample(1:n, replace=TRUE, size=n)
  xstar <- x[i]
  ystar <- y[i]
  Lb <- lm(ystar ~ xstar)
  mse = mean(Lb$residuals^2)
 })

bootCases <- out
bootCases[1:20]

# 9.6 Refer to the mammals data in the MASS package, discussed in this chapter. 
# Use your solution to the previous problem to bootstrap the MSE of the model (9.4). 
# Using the jackknife-after-bootstrap, identify which points are influential. 
# Compare this with the influential points identified from the bootstrapped slopes.
rm(list = ls());set.seed(1)
m <- 2000
data(mammals,package = 'MASS')
y <- log(mammals$brain)
x <- log(mammals$body)
n <- NROW(x)
out <- replicate(m, expr={
  i <- sample(1:n, replace=TRUE, size=n)
  xstar <- x[i]
  ystar <- y[i]
  Lb <- lm(ystar ~ xstar)
  mse = mean(Lb$residuals^2)
})

bootCases <- out
bootCases[1:20]

# identify which points are influential
rm(list = ls());set.seed(1)
library(boot);library(MASS)
theta_boot <- function(dat, ind) {
  # function to compute the patch ratio statistic
  y <- log(dat$brain)[ind]
  x <- log(dat$body)[ind]
  Lb <- lm(y ~ x)
  mse = mean(Lb$residuals^2)
}
boot.out <- boot(mammals, theta_boot, R = 2000)
jack.after.boot(boot.out, useJ=TRUE, stinf=FALSE)
jack.after.boot(boot.out, useJ=TRUE, stinf=T)

# 9.7 Plot the bootstrapped intercepts from Example 9.9 using the plot method 
# for boot with jack=TRUE. This displays a histogram of empirical influence values with a Q-Q plot, 
# and a plot similar to jack.after. boot below. Identify influential observations from the plot. 
# Are there points with standardized influence values larger than 2 in absolute value? 
# Repeat this for the slopes by setting index = 2. Do you find the same points influential in both plots?


               

