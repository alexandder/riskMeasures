library('fExtremes')

#wczytywanie danych

s1 <- unlist(read.csv("msft.csv",)[2], use.names = FALSE)
s2 <- unlist(read.csv("appl.csv")[2], use.names = FALSE)
s3 <- unlist(read.csv("yhoo.csv")[2], use.names = FALSE)
s4 <- unlist(read.csv("ibm.csv")[2], use.names = FALSE)

# initial capital
V <- 1000000
weight <- 0.25


totalNumber <- length(s1) 
numberToNovember <- length(s1) - 22
numberInDecember <- 22

alpha <- 0.95

# number of shares
lambda1 <- round(V*weight/s1[totalNumber])
lambda2 <- round(V*weight/s2[totalNumber])
lambda3 <- round(V*weight/s3[totalNumber])
lambda4 <- round(V*weight/s4[totalNumber])

lambda <- c(lambda1, lambda2, lambda3, lambda4)

calculateLogReturns <- function(prices) {
  return(diff(log(prices), lag=1))
}

# matrix of prices
S <- cbind(s1, s2, s3, s4)

# matrix of log returns
X <- cbind(calculateLogReturns(s1), calculateLogReturns(s2), calculateLogReturns(s3), calculateLogReturns(s4))


expMinusOne <- function(x) {
  return(exp(x) - 1)
}

# nonlinear loses
L <- c()

for (i in 1:totalNumber-1) {
  L[i] <- -(lambda*S[i,]) %*% expMinusOne(X[i,])
}

sre = function(n=totalNumber-1, f=L){
  x <- sort(L)
  u <- array(0, c(n-1))
  for (i in 1:(n-1)) {
    u[i] <- mean(x[(i+1):n]) - x[i]
  }
  return(list(x=x[2:n], u=u))
}

# Rysowanie wykresu (xi,ei)

p=sre()
plot(p$x, p$u)

u <- 5000

# dopasowanie uog?lnionego rozk?adu Pareto do danych d

fit = gpdFit(L, u, type = "pwm")
print(fit)

xi <- -0.116
beta <- 5314.417

numberOverThreshold <- length(subset(L, L>u))

VAR95 <- beta/xi * (((1-alpha)*(totalNumber/numberOverThreshold))^(-xi) - 1) + u
AVAR95 <- VAR95/(1-xi) + (beta-xi*u)/(1-xi)


#backtesting

LosesToNovember <- c()

for (i in 1 : numberToNovember-1) {
  LosesToNovember[i] <- -(lambda*S[i,]) %*% expMinusOne(X[i,])
}

p2 <- sre(numberToNovember-1, LosesToNovember)
plot(p2$x, p2$u)

fit2 = gpdFit(LosesToNovember, u, type = "pwm")
print(fit2)

xi2 <- -0.133
beta2 <- 5490.218

numberOverThresholdNovember <- length(subset(LosesToNovember, LosesToNovember>u))

VARToNovember <- beta2/xi2 * (((1-alpha)*(numberToNovember/numberOverThresholdNovember))^(-xi2) - 1) + u

LosesInDecember <- c()
for (j in 1 : numberInDecember-1) {
  LosesInDecember[j] <- -(lambda*S[numberToNovember+j,]) %*% expMinusOne(X[numberToNovember+j,])
}

#wektor 0-1; 1 to znaczy, ze jakas strata w grudniu byla wieksza niz VAR z calego okresu
It <- sapply(LosesInDecember>VARToNovember, as.integer)

binom.test(sum(It), length(It), 1-alpha, alternative = "two.sided")