### Sheet 2 ###
#Ex. 1
#From the lecture, we know that the posterior is distributed as beta(alpha_0+k,beta_0+n-k) if the prior is beta(alpha_0,beta_0)
#Consequently, we simulate n draws from a Bernoulli distribution with fixed parameter theta. Then, the number of draws n 
#is increased to investigate the asymptotic properties of the posterior.
library(MASS)
theta <- 0.3
alpha_0 <- 1
beta_0  <- 1  ##this gives an uninformative (flat) prior

n <- 10
k <- sum(rbinom(n = 10, size = 1, prob = theta))

# draw a sample from the posterior distribution of theta
post_sample <- rbeta(1000, alpha_0 + k, beta_0 + n - k)
truehist(post_sample)


#now increase n in a sequence from 10 to 100000
n_seq    <- c(10,20,40,100,200,500,1000,5000,10000,100000)
normgrid <- seq(from=-4, to=4, length=1000)
for(i in 1:length(n_seq)){
  k <- sum(rbinom(n = n_seq[i], size = 1, prob = theta))                #draw a sample form the bernoulli distribution
  post_sample       <- rbeta(1000, alpha_0 + k, beta_0 + n_seq[i] - k)  #draw a sample from the posterior
  truehist(post_sample, main=paste("Sample size:", n_seq[i]), xlab="Posterior: p(theta|y)") #draw a histogram of the posterior
  ML_est            <- k/n_seq[i]                                       #compute the ML estimate
  abline(v=ML_est, col="red", lwd=2)                                    #add a vertical line at the ML estimate to the histogram
  abline(v=theta, col="green", lwd=2)                                   #add a vertical line at the value used for simulation to the histogram
  standardized_post <- (post_sample-mean(post_sample))/sd(post_sample)  #compute the standardized posterior
  readline()                                                            #wait for user input (ie "enter")    
  truehist(standardized_post, main=paste("Standardized Posterior, sample size:", n_seq[i])) #draw the histogram of the standardied posterior
  lines(normgrid, dnorm(normgrid), col="red")                           #draw the density of a std. normal distribution into the histogram
  readline()                                                            #wait for user input (ie "enter")    
}




###Sheet 4
#Ex 3
#compute the quantiles of the exponential distribution
my_qexp <- function(p, lambda = 3){
  return(log(1-p)/(-lambda))
}
#create a sample from the exponential distribution with rate = 3
exp_3_sample <- my_qexp(runif(5000))

#Ex 4
#define the non normalized density, ie it integrates not to one
prop_dens <- function(x){  
  return(cos(exp(-(x-2)^2)))*(x<=4)*(x>=0)
}

#compute the proportionality constant
prop_const <- integrate(prop_dens, 0,4)

#define the density
dens <- function(x){
  return(prop_dens(x)/3.409319)
}

#define the cdf
prob <- function(x){
  return(integrate(dens, lower=0, upper=x)$value)
}

#define a function to invert the cdf using optim
to_min <- function(x, q){
  return(abs(q-prob(x)))
}

#invert the cdf
quant <- function(q){
  return(optimise(to_min, interval = c(0,4), q=q)$minimum)
}

#draw a uniformly distributed random sample and initialize an empty vector to store the sample in
quantiles <- runif(5000)
sampled    <- rep(0,5000)
#use the inversion method to sample from the given distribution
for(i in 1:5000){
  sampled[i] <- quant(quantiles[i])
}
truehist(sampled)


#Ex 6
# The given distribution function is the one of the Clayton Copula with parameter theta=0.5.
#To sample from this distribution we first sample a value for X, then look at the conditional distribution of Y,
# given the sampled value of X and draw a sample from this conditional distribution using the standard inversion technique.

my_Clayton_cond_inverse <- function(u,p){
  return(cbind(u,((p^(1/3)*sqrt(u))/(1-p^(1/3)+p^(1/3)*sqrt(u)))^2))
}
N <- 10000
clayton_sample <- my_Clayton_cond_inverse(runif(N), runif(N))

cov(clayton_sample)
sum((clayton_sample[,1]>0.4&clayton_sample[,1]<0.6&clayton_sample[,2]>0.4&clayton_sample[,2]<0.6))/N
sum(apply(clayton_sample,1,sum)<=0.5)/N

#Ex 7
# First we need to compute the maximum of the given "density". Using the maximum we can create a rectangle that just covers
# the density at it's maximum. Then we sample points from the rectangle and reject those that are above the density.
# The maximum can be found analytically, it is located at 1/5 and has a value of approx. 0.082. Hence, a rectangle between
# the points (0, 0), (1, 0), (0, 0.083), (1, 0.083) is suited to sample from the given density.

beta_no_dens <- function(x,alpha,beta){
  return(x^(alpha-1)*(1-x)^(beta-1))
}

x_cords     <- runif(100000)
y_cords     <- runif(100000, max = 0.083)
sample      <- x_cords[beta_no_dens(x_cords,2,5)>y_cords]
truehist(sample)
length(sample)


mean(sample)    #theoretical value: alpha/(alpha+beta)=2/(2+5)=2/7~0.2857 -> very accurate
var(sample)     #theoretical value: (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))~0.02551 -> very accurate



###GitKrakentestcomment

#just another comment
