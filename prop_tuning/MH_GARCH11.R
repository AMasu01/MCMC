GARCH_MH <- function(data, adaption = 1, steps = 2000, burn_in = 1000, tuning = 1000, proposal = 1){
  steps      <- steps + burn_in + tuning
  MH_mat     <- matrix(nrow=steps, ncol=3)
  MH_mat[1,] <- c(a0 = 0, a1 = 0.05, b1 = 0.9)
  lnL        <- rep(0,steps)
  lnL[1]     <- loglikelihood_GARCH11(data = data, theta = MH_mat[1,])
  
  if(proposal == 1){   ###large step independent proposal
    for(i in 2:steps){
      theta_proposal <- MH_mat[i-1,]+rnorm(3, sd = 0.001)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = theta_proposal)
      R <- log(runif(1))
      if((lnL_proposal-lnL[i-1])>R){
        MH_mat[i,] <- theta_proposal
        lnL[i]     <- lnL_proposal
      }else{
        MH_mat[i,] <- MH_mat[i-1,]
        lnL[i]     <- lnL[i-1]
      }
    }
  }
  
  
  if(proposal == 2){   ###small step independent proposal
    for(i in 2:steps){
      theta_proposal <- MH_mat[i-1,]+rnorm(3, sd = 0.000001)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = theta_proposal)
      R <- log(runif(1))
      if((lnL_proposal-lnL[i-1])>R){
        MH_mat[i,] <- theta_proposal
        lnL[i]     <- lnL_proposal
      }else{
        MH_mat[i,] <- MH_mat[i-1,]
        lnL[i]     <- lnL[i-1]
      }
    }
  }
  
  
  if(proposal == 3){   ###individual step independent proposal
    for(i in 2:steps){
      theta_proposal <- MH_mat[i-1,]+rnorm(3, sd = c(1e-6,0.01,0.01))
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = theta_proposal)
      R <- log(runif(1))
      if((lnL_proposal-lnL[i-1])>R){
        MH_mat[i,] <- theta_proposal
        lnL[i]     <- lnL_proposal
      }else{
        MH_mat[i,] <- MH_mat[i-1,]
        lnL[i]     <- lnL[i-1]
      }
    }
  }
  
  
  
  if(proposal == 4){    ### three sampling blocks! 
    for(i in 2:steps){
      a0_proposal    <- MH_mat[i-1,1]+rnorm(1, sd=1e-6)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(a0_proposal, MH_mat[i-1,2:3]))
      R <- log(runif(1))
      if((lnL_proposal-lnL[i-1])>R){
        MH_mat[i,] <- c(a0_proposal, MH_mat[i-1,2:3])
        lnL[i]     <- lnL_proposal
      }else{
        MH_mat[i,] <- MH_mat[i-1,]
        lnL[i]     <- lnL[i-1]
      }
      
      a1_proposal    <- MH_mat[i-1,2] + rnorm(1, sd=0.01)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(MH_mat[i,1], a1_proposal, MH_mat[i-1,3]))
      R <- log(runif(1))
      if((lnL_proposal-lnL[i])>R){
        MH_mat[i,2]  <- a1_proposal
        lnL[i]       <- lnL_proposal
      }
      
      b1_proposal    <- MH_mat[i-1,3] + rnorm(1, sd=0.01)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(MH_mat[i,1:2], b1_proposal))
      R <- log(runif(1))
      if((lnL_proposal-lnL[i])>R){
        MH_mat[i,3]  <- b1_proposal
        lnL[i]       <- lnL_proposal
      }
    }
  }
  
  
  if(proposal == 5){    ### three sampling blocks with tuning of the proposal (just variances)
    #tuning (and also some burn in)
    ARmat <- matrix(0, nrow=steps, ncol=3)
    sds   <- c(1e-6, 0.01, 0.01)
    for( j in 2:tuning){
      if(j>101){
        ARmat[j,1] <- sum(diff(MH_mat[(j-100):(j-1),1])!=0)/99
        ARmat[j,2] <- sum(diff(MH_mat[(j-100):(j-1),2])!=0)/99
        ARmat[j,3] <- sum(diff(MH_mat[(j-100):(j-1),3])!=0)/99
        if(j%%50==0){
          if(mean(ARmat[(j-49):j,1])<0.2){sds[1] <- sds[1]*0.9}else{if(mean(ARmat[(j-49):j,1])>0.3){sds[1] <- sds[1]*1.1}}
          if(mean(ARmat[(j-49):j,2])<0.2){sds[2] <- sds[2]*0.9}else{if(mean(ARmat[(j-49):j,2])>0.3){sds[2] <- sds[2]*1.1}}
          if(mean(ARmat[(j-49):j,3])<0.2){sds[3] <- sds[3]*0.9}else{if(mean(ARmat[(j-49):j,3])>0.3){sds[3] <- sds[3]*1.1}}
        }
      }
      a0_proposal    <- MH_mat[j-1,1]+rnorm(1, sd = sds[1])
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(a0_proposal, MH_mat[j-1,2:3]))
      R              <- log(runif(1))
      if((lnL_proposal-lnL[j-1])>R){
        MH_mat[j,] <- c(a0_proposal, MH_mat[j-1,2:3])
        lnL[j]     <- lnL_proposal
      }else{
        MH_mat[j,] <- MH_mat[j-1,]
        lnL[j]     <- lnL[j-1]
      }
      
      a1_proposal    <- MH_mat[j-1,2] + rnorm(1, sd = sds[2])
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(MH_mat[j,1], a1_proposal, MH_mat[j-1,3]))
      R <- log(runif(1))
      if((lnL_proposal-lnL[j])>R){
        MH_mat[j,2]  <- a1_proposal
        lnL[j]       <- lnL_proposal
      }
      
      b1_proposal    <- MH_mat[j-1,3] + rnorm(1, sd = sds[3])
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(MH_mat[j,1:2], b1_proposal))
      R <- log(runif(1))
      if((lnL_proposal-lnL[j])>R){
        MH_mat[j,3]  <- b1_proposal
        lnL[j]       <- lnL_proposal
      }
    }
    
    
    #Burn-in and sampling 
    for(i in (j+1):steps){
      if(i>101){
        ARmat[i,1] <- sum(diff(MH_mat[(i-100):(i-1),1])!=0)/99
        ARmat[i,2] <- sum(diff(MH_mat[(i-100):(i-1),2])!=0)/99
        ARmat[i,3] <- sum(diff(MH_mat[(i-100):(i-1),3])!=0)/99
      }
      a0_proposal    <- MH_mat[i-1,1]+rnorm(1, sd = sds[1])
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(a0_proposal, MH_mat[i-1,2:3]))
      R <- log(runif(1))
      if((lnL_proposal-lnL[i-1])>R){
        MH_mat[i,] <- c(a0_proposal, MH_mat[i-1,2:3])
        lnL[i]     <- lnL_proposal
      }else{
        MH_mat[i,] <- MH_mat[i-1,]
        lnL[i]     <- lnL[i-1]
      }
      
      a1_proposal    <- MH_mat[i-1,2] + rnorm(1, sd = sds[2])
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(MH_mat[i,1], a1_proposal, MH_mat[i-1,3]))
      R <- log(runif(1))
      if((lnL_proposal-lnL[i])>R){
        MH_mat[i,2]  <- a1_proposal
        lnL[i]       <- lnL_proposal
      }
      
      b1_proposal    <- MH_mat[i-1,3] + rnorm(1, sd = sds[3])
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = c(MH_mat[i,1:2], b1_proposal))
      R <- log(runif(1))
      if((lnL_proposal-lnL[i])>R){
        MH_mat[i,3]  <- b1_proposal
        lnL[i]       <- lnL_proposal
      }
    }
  }
  
  
  
  if(proposal == 6){    ### single sampling block with tuned covariance matrix
    if(tuning < 2000){tuning <- 2000}  #at least 2000 tuning steps 
    library(copula)
    library(mvtnorm)
    #tuning (and also some burn in)
    ARatio    <- rep(0, steps)
    sds       <- c(1e-6, 0.01, 0.01)
    rho       <- rep(0,3)
    for( j in 2:tuning){
      if(j>101){
        ARatio[j] <- sum(diff(MH_mat[(j-100):(j-1),1])!=0)/99
        if(j%%50==0){
          if(mean(ARatio[(j-49):j])<0.3){sds <- sds * 0.9}else{if(mean(ARatio[(j-49):j])>0.4){sds <- sds * 1.1}}
        }
      }
      if(i>=(ceiling(tuning/2)-1)&i%%500==0){
        rho <- P2p(cor(MH_mat[(j-500):(j-1)]))
      }
      covmat         <- (sds%*%t(sds))*p2P(rho)
      theta_proposal <- MH_mat[j-1,] + rmvnorm(1, sigma = covmat)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = theta_proposal)
      R              <- log(runif(1))
      if((lnL_proposal-lnL[j-1])>R){
        MH_mat[j,] <- theta_proposal
        lnL[j]     <- lnL_proposal
      }else{
        MH_mat[j,] <- MH_mat[j-1,]
        lnL[j]     <- lnL[j-1]
      }
    }
    
    
    #Burn-in and sampling 
    for(i in (j+1):steps){
      if(i>101){
        ARatio[i] <- sum(diff(MH_mat[(i-100):(i-1),1])!=0)/99
      }
      theta_proposal <- MH_mat[i-1,] + rmvnorm(1, sigma = covmat)
      lnL_proposal   <- loglikelihood_GARCH11(data = data, theta = theta_proposal)
      R              <- log(runif(1))
      if((lnL_proposal-lnL[i-1])>R){
        MH_mat[i,] <- theta_proposal
        lnL[i]     <- lnL_proposal
      }else{
        MH_mat[i,] <- MH_mat[i-1,]
        lnL[i]     <- lnL[i-1]
      }
    }
    return(list(MH_mat[(1+burn_in+tuning):N,], lnL[(1+burn_in+tuning):N], ARatio[(1+burn_in+tuning):N]))  
    
  }
  
  
  
  #return(list(MH_mat, lnL, ARmat))  
  
  return(list(MH_mat[(1+burn_in+tuning):N,], lnL[(1+burn_in+tuning):N], ARmat[(1+burn_in+tuning):N,]))  
  
}