#### Proposal tuning ####
#A GARCH(1,1)-model is to be estimated to illustrate the importance ofan adequate proposal distribution.
# x_t = gamma_0 + sigma_t *eps_t; eps_t~N(0,1)
# sigma_t^2 = a_0 + a_1*x_(t-1) + b_1 *sigma_(t-1)^2

set.seed(1234)
theta    <- c(gamma_0 = 0, a0 = 1e-6, a1 = 0.07, b1 = 0.92)
N        <- 1000
x        <- rep(0,N)
x[1]     <- 0.1
sigma2   <- rep(0,N)
sigma2[1]<- 1e-4
plot(x, type="l", main="Daily return", xlab="day", ylab = "return")
plot(sigma2, type="l", main="Volatility", xlab="day", ylab="volatilty")


for(i in 2:N){
  sigma2[i] <- theta[2] + theta[3]*x[i-1]^2 + theta[4]*sigma2[i-1]
  x[i]      <- theta[1] + sqrt(sigma2[i])*rnorm(1) 
}






