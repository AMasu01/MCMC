loglikelihood_GARCH11 <- function(data, theta){
  gamma_0   <- 0
  a0        <- theta[1]
  a1        <- theta[2]
  b1        <- theta[3]
  N         <- length(data)
  sigma2    <- rep(0,N)
  sigma2[1] <- a0/(1-a1-b1)
  eps       <- rep(1,N)
  if((abs(a1)+abs(b1))>=1){
    #print("bla")
    return(-666666)
  }
  if(any(theta[-1]<0)){return(-666666)}
  for(i in 2:N){
    sigma2[i] <- a0 + a1 * data[i-1]^2 + b1 * sigma2[i-1]
    if(sigma2[i]<0){return(-666666)}
    eps[i]    <- data[i]/sqrt(sigma2[i])
  }
  #plot(sigma2, type = "l")
  #plot(sigma2, type="l")
  #plot(sigma2[500:1000], type="l")
  #plot(eps, type="l")
  #truehist(eps)
  lnL <- sum(dnorm(eps[2:N], log = T))-sum(log(sigma2[2:N]))
  if(is.na(lnL)|is.infinite(lnL)){return(-666666)}
  return(lnL)
}
