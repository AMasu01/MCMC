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
