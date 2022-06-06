
################################################################################
################################################################################

## binomial model - with less than threshold 

## generate data: includes intercept 

gen_data_arm <- function(n=c(300),rand_pr=NULL,n_arm=NULL){
  ##
  ## n is the sample in one interim 
  ##
  dat <- data.table::CJ(
    id = 1:n,
    id_arm = NA_integer_
  )
  ##
  if(is.null(rand_pr)){
    stop("provide value for rand_pr")
  }
  ##
  if(is.null(n_arm)){
    n_arm <- length(rand_pr)
    dat$id_arm <- sample(1:n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  }
  else{
    dat$id_arm <- sample(n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  }
  ##  
  dat
}
#input <- gen_data_arm()

##
#c(0,-1,0.5,0) # beta
#p = 0.5 # values below threshold
#b0 <- qlogis(p) #
#p1 <- plogis(b0 - 1)
#p2 <- plogis(b0 + 0.5)
#p3 <- plogis(b0 + 0)
#rbinom(25,1,prob=p)
#rbinom(25,1,prob=p1)
#rbinom(25,1,prob=p2)
#rbinom(25,1,prob=p3)
##
gen_data_bino_copd <- function(N, prob_below_threshold = NULL,
                               bta=c(0,-1,0.5,0), 
                               rand_pr=c(0.25,0.25,0.25,0.25)){
  ##
  input <- data.frame(gen_data_arm(n=N, rand_pr = rand_pr, n_arm = NULL))
  n <- table(input$id_arm)
  ##
  K <- length(bta)
  p <- c()
  input$y <- NA
  for(k in as.numeric(names(n))){
    if(k==1){
      if(is.null(prob_below_threshold)){
        p[k] <- plogis(bta[k])
      }
      else{
        bta[k] <- qlogis(prob_below_threshold)
        p[k] <- prob_below_threshold
      }
      input[input$id_arm==k,"y"] <- rbinom(n[names(n)%in%k],1,prob=p[k])
    }
    else{
      p[k] <- plogis(bta[1]+bta[k])
      input[input$id_arm==k,"y"] <- rbinom(n[names(n)%in%k],1,prob=p[k])
    }
  }
  #print(paste0("Pr(obs below the threshold) = ", mean(input$y)))
  ##
  out <- data.frame(input, x=as.factor(input$id_arm))
  out
}
##
dat <- gen_data_bino_copd(N=100, prob_below_threshold = 0.2)
##

##
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)  
##
model_bernoulli <- stan_model(file="ber_logit.stan")
##
trial_fast_bernoulli <- function(data=gen_data_bino_copd(N=100, 0.2), 
                       model=model_bernoulli,
                       mu_beta = c(0,0,0,0), # with intercept -
                       sig_beta = c(5,5,5,5)) # with intercept - 
{
  ##
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)  
  X <- model.matrix( ~ as.factor(x), data = data)
  y <- data$y
  K <- ncol(X) # with interpect
  N <- nrow(X)
  if(!length(mu_beta)%in%K){
    stop("check mu_beta")
  }
  if(!length(sig_beta)%in%K){
    stop("check sig_beta")
  }
  ##
  berdat <- list(
      N = N,
      K = K, # with interpect
      x = X,
      y = y,
      mu_beta = mu_beta, sig_beta=sig_beta
    )
  ##
  ber.fast <- optimizing(model, data=berdat, hessian = TRUE)
  # Covariance is negative inverse Hessian
  opt_cov <- solve(-ber.fast$hessian)
  ##
  post.beta.mn <- ber.fast$par[1:K]
  post.beta.sd <- sqrt(diag(opt_cov[1:K,1:K]))
  ##
  rm(ber.fast)
  ##
  sample_size <- unique(sapply(data,length))
  ##
  return(list(post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              sample_size=sample_size))
}
##
## trail interim - COPD
##
trail_interim_bino_copd <- function(model, 
                              prob_below_threshold = NULL, # i.e.,
                              n_sample=c(100,50,50), # n_sample => increment values
                              rand_assignment=c(1,1,1,1),   
                              true_beta = c(0,-1,0.5,0), # first input is intercept
                              mu_beta = c(0,0,0,0), # with intercept -
                              sig_beta = sqrt(c(5,5,5,5)), # with intercept - 
                              D0=0.10, D1=0.90,
                              capType="noCap"){
  ##
  ##
  K <- length(true_beta)
  n_j <- length(n_sample)
  y_init <- list()
  n_store <- c()
  per_below_thrs <- matrix(NA,K,n_j) # arms specific missing
  n_interimSS <- matrix(NA,K,n_j)
  n_rand_pr <- matrix(NA,K,n_j)
  pr.decision <- matrix(NA,K-1,n_j)
  out <- list()
  for(j in 1:n_j){
    if(j==1){
      N <- n_sample[j]
      rand_pr <- rand_assignment/sum(rand_assignment)
      y_init[[j]] = gen_data_bino_copd(N=N, prob_below_threshold = prob_below_threshold, 
                                       bta=true_beta, rand_pr=rand_pr)[,3:4] # only taking the y and x variables 
      n_store[j] <- nrow(y_init[[j]])
      n_interimSS[,j] <- table(y_init[[j]]$x)
      n_rand_pr[,j] <- rand_pr
      #per_below_thrs[j] <- mean(y_init[[j]]$y)
      per_below_thrs[,j] <- table(y_init[[j]][y_init[[j]]$y%in%1,])/n_interimSS[,j] # armwise missing 
      ##
      out[[j]] = trial_fast_bernoulli(data = y_init[[j]], 
                            model = model, # stan model from stan_model code 
                            mu_beta = mu_beta, # with intercept -
                            sig_beta = sig_beta)
      ck <- NULL
      for(l in 2:K){
        #ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))<1))
        ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))>1))
      }
      pr.decision[,j] <- ck
      ##
    }
    else{
      N <- n_sample[j]
      ## For no sample size cap: if one arm is dropped then other arms 
      ## may receive more participants up to the nMax[sum(n_sample)] participants 
      if(capType=="noCap"){
        drop_arm <- c(which(pr.decision[,j-1] < D0) + 1,which(pr.decision[,j-1] > D1) + 1)
        rand_pr[drop_arm] <- 0
        rand_pr <- rand_pr/sum(rand_pr)
        y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],
                                        gen_data_bino_copd(N=N, prob_below_threshold = prob_below_threshold, bta=true_beta, rand_pr=rand_pr)[,3:4],
                                        SIMPLIFY=FALSE))
        n_store[j] <- nrow(y_init[[j]])
        n_interimSS[,j] <- table(y_init[[j]]$x)
        n_rand_pr[,j] <- rand_pr
        #per_below_thrs[j] <- mean(y_init[[j]]$y)
        per_below_thrs[,j] <- table(y_init[[j]][y_init[[j]]$y%in%1,])/n_interimSS[,j] # armwise missing 
      }
      else{
        y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],
                                        gen_data_bino_copd(N=N, prob_below_threshold = prob_below_threshold, bta=true_beta, rand_pr=rand_pr)[,3:4],
                                        SIMPLIFY=FALSE))
        n_store[j] <- nrow(y_init[[j]])
        n_interimSS[,j] <- table(y_init[[j]]$x)
        n_rand_pr[,j] <- rand_pr
        #per_below_thrs[j] <- mean(y_init[[j]]$y)
        per_below_thrs[,j] <- table(y_init[[j]][y_init[[j]]$y%in%1,])/n_interimSS[,j] # armwise missing 
      }
      ##
      out[[j]] = trial_fast_bernoulli(data = y_init[[j]], 
                                      model = model, # stan model from stan_model code 
                                      mu_beta = mu_beta, # with intercept -
                                      sig_beta = sig_beta)
      ck <- NULL
      for(l in 2:K){
        #ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))<1))
        ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))>1))
      }
      pr.decision[,j] <- ck
      ##
    }
  }
  ##
  post.beta.mn <- sapply(out, function(x) x$post.beta.mn)
  post.beta.sd <- sapply(out, function(x) x$post.beta.sd)
  ##
  dimnames(post.beta.mn)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(post.beta.sd)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(pr.decision)[[1]] = c(paste0("treat",1:(length(true_beta)-1)))
  ##
  dimnames(n_interimSS)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_interimSS)[[2]] = c(paste0("interim",1:n_j))
  ##
  dimnames(n_rand_pr)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_rand_pr)[[2]] = c(paste0("interim",1:n_j))
  ##
  per_below_thrs <- per_below_thrs*100
  dimnames(per_below_thrs)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(per_below_thrs)[[2]] = c(paste0("interim",1:n_j))
  #
  return(list(pr.decision=pr.decision,
              post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              sample_size=n_interimSS, rand_pr=n_rand_pr, 
              n_store=n_store, per_below_thrs=per_below_thrs))
}
##
## trial simulator - COPD
##
trial_interim_simulator_bino_copd <- function(nSim=10,model, 
                                              prob_below_threshold = NULL, # i.e.,
                                              n_sample=c(100,50,50), # n_sample => increment values
                                              rand_assignment=c(1,1,1,1),   
                                              true_beta = c(0,-1,0.5,0), # first input is intercept
                                              prior_mu_beta = c(0,0,0,0), # with intercept - hyper para
                                              prior_sig_beta = sqrt(c(5,5,5,5)), # with intercept - hyper para
                                              D0=0.10, D1=0.90,
                                              capType="noCap"){
  ##
  library(pbapply)
  out <- pblapply(1:nSim, 
                  function(x) trail_interim_bino_copd(model = model, 
                                                      prob_below_threshold = prob_below_threshold, # i.e.,
                                                      n_sample = n_sample, 
                                                      rand_assignment = rand_assignment,   
                                                      true_beta = true_beta, 
                                                      mu_beta = prior_mu_beta, 
                                                      sig_beta = prior_sig_beta, 
                                                      D0 = D0, D1 = D1,
                                                      capType = capType))
  ##
  n_store <- sapply(out, function(x) x$n_store)
  if(is.matrix(n_store)){
    dimnames(n_store)[[1]] = paste0("interim",1:length(n_sample))
  }
  else{
    n_store <- as.matrix(n_store)
    n_store <- t(n_store)
    dimnames(n_store)[[1]] = paste0("interim",1:length(n_sample))
  }
  ##
  tmp <- sapply(out, function(x) x$pr.decision)
  pr.decision <- array(tmp,dim=c(length(prior_mu_beta)-1,length(n_sample),nSim))
  dimnames(pr.decision)[[1]] = paste0("treat",1:(length(prior_mu_beta)-1))
  dimnames(pr.decision)[[2]] = paste0("interim",1:length(n_sample))
  ##
  tmp <- sapply(out, function(x) x$post.beta.mn)
  post.beta.mn <- array(tmp,dim=c(length(prior_mu_beta),length(n_sample),nSim))
  dimnames(post.beta.mn)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
  dimnames(post.beta.mn)[[2]] = paste0("interim",1:length(n_sample))
  ##
  tmp <- sapply(out, function(x) x$post.beta.sd)
  post.beta.sd <- array(tmp,dim=c(length(prior_mu_beta),length(n_sample),nSim))
  dimnames(post.beta.sd)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
  dimnames(post.beta.sd)[[2]] = paste0("interim",1:length(n_sample))
  ##
  sample_size <- sapply(out, function(x) x$sample_size)
  sample_size <- array(sample_size,dim=c(length(prior_mu_beta),length(n_sample),nSim))
  dimnames(sample_size)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
  dimnames(sample_size)[[2]] = paste0("interim",1:length(n_sample))
  ##
  rand_pr <- sapply(out, function(x) x$rand_pr)
  rand_pr <- array(rand_pr,dim=c(length(prior_mu_beta),length(n_sample),nSim))
  dimnames(rand_pr)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
  dimnames(rand_pr)[[2]] = paste0("interim",1:length(n_sample))
  ##
  per_below_thrs <- sapply(out, function(x) x$per_below_thrs)
  per_below_thrs <- array(per_below_thrs,dim=c(length(prior_mu_beta),length(n_sample),nSim))
  dimnames(per_below_thrs)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
  dimnames(per_below_thrs)[[2]] = paste0("interim",1:length(n_sample))
  #per_below_thrs <- sapply(out, function(x) x$per_below_thrs)
  #if(is.matrix(per_below_thrs)){
  #  dimnames(per_below_thrs)[[1]] = paste0("interim",1:length(n_sample))
  #}
  #else{
  #  per_below_thrs <- as.matrix(per_below_thrs)
  #  per_below_thrs <- t(per_below_thrs)
  #  dimnames(per_below_thrs)[[1]] = paste0("interim",1:length(n_sample))
  #}
  ##
  return(list(sample_size=sample_size, 
              rand_pr=rand_pr,
              n_store=t(n_store),
              per_below_thrs=per_below_thrs,
              pr.decision=pr.decision,
              post.beta.mn=post.beta.mn, 
              post.beta.sd=post.beta.sd,
              n_sample=n_sample, nSim=nSim,
              true_beta=true_beta, 
              prior_mu_beta=prior_mu_beta, prior_sig_beta=prior_sig_beta,
              capType=capType, D0=D0, D1=D1))
  ##
}
##
##

## simulation (#1)

model_bernoulli <- stan_model(file="ber_logit.stan")

#true_bta = cbind(rep(0,6),c(0,-0.5,-1,-1.5,-2,-2.5),c(0,-0.5,-1,-1.5,-2,-2.5),c(0,-0.5,-1,-1.5,-2,-2.5))
true_bta = cbind(rep(0,6),c(0,0.25,0.5,1,1.5,2),c(0,0.25,0.5,1,1.5,2),c(0,0.25,0.5,1,1.5,2))
true_pr = seq(0.1,0.4,0.05)
nSim = 1000
n_sample = c(100,50,50)
rand_assignment = c(1,1,1,1)
D0 = 0.1; D1 = 0.90
store_res <- NULL
id_scenario <- 0
for(pr in 1:length(true_pr)){
  for(bt1 in 1){
    for(bt2 in 1:nrow(true_bta)){
      for(bt3 in 1){#1:nrow(true_bta)){
        for(bt4 in 1){#1:nrow(true_bta)){
          res <- trial_interim_simulator_bino_copd(nSim=nSim,model=model_bernoulli, 
                                                   prob_below_threshold = true_pr[pr], # i.e.,
                                                   n_sample=n_sample, # n_sample => increment values
                                                   rand_assignment=c(1,1,1,1),   
                                                   true_beta = c(true_bta[bt1,1],true_bta[bt2,2],true_bta[bt3,3],true_bta[bt4,4]), # first input is intercept
                                                   prior_mu_beta = c(0.5,0,0,0), # with intercept - hyper para
                                                   prior_sig_beta = sqrt(c(5,5,5,5)), # with intercept - hyper para
                                                   D0=D0, D1=D1,
                                                   capType="noCap")
          
          ##
          ## summary stat
          ##
          dat <- data.table::CJ(
            interim = 1:length(n_sample),
            arm = 1:length(rand_assignment),
            id_scenario = NA_integer_,
            true_pr_missing = NA_real_,
            per_missing = NA_real_,
            true_beta = NA_real_,
            est_beta = NA_real_,
            est_exp_beta = NA_real_,
            est_sample =  NA_real_,
            est_sd =  NA_real_
          )
          id_scenario <- id_scenario + 1
          dat$id_scenario <-  id_scenario
          dat$true_pr_missing <- true_pr[pr]
          dat$true_beta = rep(c(qlogis(true_pr[pr]),true_bta[bt2,2],true_bta[bt3,3],true_bta[bt4,4]),length(n_sample))
          for(j in 1:length(n_sample)){
            #
            dat[dat$interim==j,"est_beta"] <- rowMeans(res$post.beta.mn[,j,])
            dat[dat$interim==j,"est_exp_beta"] <- rowMeans(exp(res$post.beta.mn[,j,]))
            #
            dat[dat$interim==j,"est_sample"] <- rowMeans(res$sample_size[,j,])
            dat[dat$interim==j,"est_sd"] <- apply(res$sample_size[,j,],1,sd)
            #
            dat[dat$interim==j,"per_missing"] <- rowMeans(res$per_below_thrs[,j,])
          }
          pr_decision <- NULL
          for(j in 1:length(n_sample)){
            pr_decision <- rbind(pr_decision,cbind(j,t(rbind(rowMeans(res$pr.decision[,j,]<D0),rowMeans(res$pr.decision[,j,]>D1)))))
          }
          dimnames(pr_decision)[[2]] <- c("interim","D0_0.1","D1_0.9")
          pr_decision <- data.frame(pr_decision)
          pr_decision$arm <- 2:4
          dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
          ##
          store_res <- rbind(store_res,dat)
          write.csv(store_res,file=paste0("result_bino_D0.10_D1.90_",Sys.Date(),".csv"),row.names=FALSE)
        }
      }
    }
  }
}

################################################################################
################################################################################




