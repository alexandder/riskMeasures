
s1 <- unlist(read.csv("msft.csv",)[2], use.names = FALSE)
s2 <- unlist(read.csv("appl.csv")[2], use.names = FALSE)
s3 <- unlist(read.csv("yhoo.csv")[2], use.names = FALSE)
s4 <- unlist(read.csv("ibm.csv")[2], use.names = FALSE)

# initial capital
V <- 1000000

weight <- 0.25

# length
n <- length(s1)

# number of shares
lambda1 <- V*weight/s1[1]
lambda2 <- V*weight/s2[1]
lambda3 <- V*weight/s3[1]
lambda4 <- V*weight/s4[1]

lambda <- c(lambda1, lambda2, lambda3, lambda4)

calculateLogReturns <- function(prices) {
  return(diff(log(prices), lag=1))
}

# matrix of prices
S <- cbind(s1, s2, s3, s4)

# matrix of log returns
X <- cbind(calculateLogReturns(s1), calculateLogReturns(s2), calculateLogReturns(s3), calculateLogReturns(s4))

# L1 <- -(lambda1*s1[1]*(exp(X[1,1]) -1) + lambda2*s2[1]*(exp(X[1,2]) - 1) + lambda3*s3[1]*(exp(X[1,3]) - 1) + lambda4*s4[1]*(exp(X[1,4]) - 1))

expMinusOne <- function(x) {
  return(exp(x) - 1)
}

# nonlinear loses
L <- c()

for (i in 1:n-1) {
  L[i] <- -(lambda*S[i,]) %*% expMinusOne(X[i,])
}

VAR99 <- quantile(L, 0.99)
VAR95 <- quantile(L, 0.95)

# Component VAR

weightComponent <- 1/3

VAR99Component <- c()
VAR95Component <- c()

for(j in 1:4) {
  # deleting column j from prices matrix
  SComponent <- S[, -j]
  
  # new number of shares
  lambdaComponent <- c(V*weightComponent/SComponent[1,1], V*weightComponent/SComponent[1,2], V*weightComponent/SComponent[1,3])
  
  # log returns matrix without j-th column
  XComponent <- X[, -j]
  
  LComponent <- c()
  
  for (k in 1:n-1) {
    LComponent[k] <- -(lambdaComponent*SComponent[k,]) %*% expMinusOne(XComponent[k,])
  }
  
  VAR99Component[j] <- quantile(LComponent, 0.99)
  VAR95Component[j] <- quantile(LComponent, 0.95)
}

# Linear approximation of loses
LLinear <- c()

# vector of weights
W <- c(weight, weight, weight, weight)

for (i in 1:n-1) {
  LLinear[i] <- -V*(W %*% X[i,])
}

VAR99LinearLoses <- quantile(LLinear, 0.99)
VAR95LinearLoses <- quantile(LLinear, 0.95)

# Marignal VAR

# covariance matrix of log returns
SigmaMatrix <- cov(X)

# little sigma; result of multiplication is 1x1 matrix so we have to take that number explicitly
sigmaSquared <- (t(W) %*% SigmaMatrix %*% W)[1,1]

VAR99Marginal <- (VAR99LinearLoses/V) * (SigmaMatrix %*% W)/sigmaSquared
VAR95Marginal <- (VAR95LinearLoses/V) * (SigmaMatrix %*% W)/sigmaSquared

# Variance-covariance method

mu <- colMeans(X)

XHat <- cbind(X[,1] - mu[1], X[,2] - mu[2], X[,3] - mu[3], X[,4] - mu[4])

SigmaHat <- 1/(n-1) * (t(XHat) %*% XHat)

a <- V*W

sigma <- sqrt(t(a) %*% SigmaHat %*% a)

v <- t(a) %*% mu

VAR99VC <- v + sigma*qnorm(0.99, mean = 0, sd=1)
VAR95VC <- v + sigma*qnorm(0.95, mean = 0, sd=1)