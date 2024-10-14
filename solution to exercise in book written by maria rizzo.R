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
