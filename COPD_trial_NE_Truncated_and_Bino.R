
################################################################################
################################################################################

library(pbapply)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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

gen_data_copd <- function(N, bta, sig, rand_pr=c(0.25,0.25,0.25,0.25)){
  ##
  input <- gen_data_arm(n=N, rand_pr = rand_pr, n_arm = NULL)
  ##
  xx <- model.matrix(~factor(input$id_arm, levels=c(1:length(rand_pr))))
  y <- xx%*%c(bta) + rnorm(N,0,sig)
  #enrol_time <- runif(N, min = 0, max = 28)
  #obs_time <- enrol_time + 7
  #out <- data.frame(enrol_time = enrol_time, obs_time = obs_time, y = y, x=as.factor(input$id_arm))
  #out[order(out$enrol_time),]
  out <- data.frame(input, y = y, x=as.factor(input$id_arm))
  out
}
#gen_data_copd(N=50,bta=c(0.5,0,0,0),sig=1,rand_pr=c(0.25,.25,.25,.25))
#gen_data_copd(N=50,bta=c(0.5,0),sig=1,rand_pr=c(0.5,0.5), n_arm=c(1,3))

##
## not using QR-decomposition
##
model_cen <- stan_model(file="lm_beta_leftcensored_optiz_knownL.stan")
model_noncen <- stan_model(file="lm_beta.stan")
##
trial_fast <- function(data, 
                       model,
                       modelType=NULL,
                       L = NULL,
                       mu_beta = c(0,0,0,0), # with intercept -
                       sig_beta = c(5,5,5,5), # with intercept - 
                       a=0.5, b=2){
  ##
  ## modelType = "nonCensored" or "Censored"
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
  if(modelType=="nonCensored"){
    lmdat <- list(
      N = N,
      K = K, # with interpect
      x = X,
      y = y,
      mu_beta = mu_beta, sig_beta=sig_beta,
      a = a, b = b
    )
  }
  else if(modelType=="Censored"){
    ##
    ## following code is for when the value is known see lm_beta_leftcensored_optiz_knownL.stan
    if (is.null(L)) {L = min(y)}
    if (min(y) < L) {stop("Minimum value below censored value")}
    if (min(y) > L) {warning("No censored values")}
    X_cens = subset(X, y == L)
    X_obs = subset(X, y > L)
    y_obs = subset(y, y > L)
    lmdat <- list(
      N_obs = nrow(X_obs),
      N_cens = nrow(X_cens),
      K = K,
      X_obs = X_obs,
      X_cens = X_cens,
      y_obs = y_obs,
      L = L, 
      mu_beta = mu_beta, sig_beta=sig_beta,
      a = a, b = b      
    )
  }
  ##
  else{
    stop("provide modeType argument")
  }
  ##
  lm.fast <- optimizing(model, data=lmdat, hessian = TRUE)
  # Covariance is negative inverse Hessian
  opt_cov <- solve(-lm.fast$hessian)
  ##
  post.beta.mn <- lm.fast$par[1:K]
  post.beta.sd <- sqrt(diag(opt_cov[1:K,1:K]))
  post.sigma.mn <- lm.fast$par[K+1]
  post.sigma.sd <- sqrt(opt_cov[K+1,K+1])
  ##
  ##
  rm(lm.fast)
  ##
  sample_size <- unique(sapply(data,length))
  ##
  return(list(post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
    post.sigma.mn=post.sigma.mn, post.sigma.sd=post.sigma.sd,
    sample_size=sample_size))
}
##
model_bino <- stan_model(file="ber_logit.stan")
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

trail_interim_copd_normal_bino <- function(model_noncen=NULL, model_cen=NULL, model_bino=NULL,
                                           modelType, L=NULL, L_per=NULL, 
                                           n_sample=c(100,50,50), # n_sample => increment values
                                           rand_assignment=c(1,1,1,1),   
                                           true_beta = c(0.5,0,0,0), # first input is intercept
                                           true_sig = c(1),
                                           mu_beta = c(0,0,0,0), # with intercept -
                                           sig_beta = c(5,5,5,5), # with intercept - 
                                           a=3, b=1,
                                           capType="cap",
                                           D0=0.10, D1=0.90){
  ##
  ## L = censored point for a left censored model 
  ## L_per = percentage of left censored observations 
  ## model_cen = use a noncensored model, useful when L yields to zero missing values in the data, not useful when we use L_per 
  ##
  K <- length(true_beta)
  n_j <- length(n_sample)
  y_init <- list()
  y_init_bino <- list()
  n_store <- c()
  per_missing <- matrix(NA,K,n_j) # arms specific missing
  n_interimSS <- matrix(NA,K,n_j)
  n_rand_pr <- matrix(NA,K,n_j)
  pr.decision <- matrix(NA,K-1,n_j)
  pr.decision_bino <- matrix(NA,K-1,n_j)
  out <- list()
  out_bino <- list()
  for(j in 1:n_j){
    if(j==1){
      N <- n_sample[j]
      rand_pr <- rand_assignment/sum(rand_assignment)
      y_init[[j]] = gen_data_copd(N=N, bta=true_beta, sig=true_sig, rand_pr=rand_pr)[,3:4] # only taking the y and x variables 
      n_store[j] <- nrow(y_init[[j]])
      n_interimSS[,j] <- table(y_init[[j]]$x)
      n_rand_pr[,j] <- rand_pr
      if(modelType=="Censored"){
        if(!is.null(L_per)){ # 
          #cen_n <- round(n_store[j]*L_per)
          cen_n <- round(n_sample[j]*L_per)
          LL <- sort(y_init[[j]]$y)[cen_n]
        }
        else{
          LL <- L
        }
        #print(LL)
        #print(L)
        y_init[[j]]$y <- ifelse(y_init[[j]]$y <= LL, LL, y_init[[j]]$y) 
        #y_init[[j]]$y[y_init[[j]]$y<=LL] <- LL
        #per_missing[j] <- mean(y_init[[j]]$y<=LL) # overall missing
        per_missing[,j] <- table(y_init[[j]][y_init[[j]]$y<=LL,"x"])/n_interimSS[,j] # armwise missing 
      }
      else{
        LL = NULL
        per_missing[,j] <- 0
      }
      ##
      #print(c(j,rand_pr))
      ##
      if(unique(per_missing[,j])[length(unique(per_missing[,j]))]%in%0){
        model = model_noncen
        modelType2 = "nonCensored"
        LL = NULL
      }
      else{
        model = model_cen
        modelType2 = "Censored"
      }
      ##
      #browser()
      #print(c(j,modelType2,per_missing[j]))
      out[[j]] = trial_fast(data = y_init[[j]], 
                            model = model, # stan model from stan_model code 
                            modelType = modelType2, L=LL,
                            mu_beta = mu_beta, # with intercept -
                            sig_beta = sig_beta, # with intercept - 
                            a=a, b=b)
      pr.decision[,j] <- pnorm(0, out[[j]]$post.beta.mn[2:K], out[[j]]$post.beta.sd[2:K])
      #drop_arm <- c(which(pr.decision[,j] < D0) + 1,which(pr.decision[,j] > D1) + 1)
      ##
      ## for bino 
      y_init_bino[[j]] = y_init[[j]]
      y_init_bino[[j]][y_init[[j]]$y<=LL,"y"] = 1
      y_init_bino[[j]][y_init[[j]]$y>LL,"y"] = 0
      out_bino[[j]] = trial_fast_bernoulli(data = y_init_bino[[j]], 
                                      model = model_bino, # stan model from stan_model code 
                                      mu_beta = mu_beta, # with intercept -
                                      sig_beta = sig_beta)
      ck <- NULL
      for(l in 2:K){
        #ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))<1))
        ck <- c(ck,mean(exp(rnorm(10000,out_bino[[j]]$post.beta.mn[l],out_bino[[j]]$post.beta.sd[l]))>1))
      }
      pr.decision_bino[,j] <- ck
      ##
    }
    else if(j==n_j){
      N <- cumsum(n_sample)[j]-n_store[j-1]
      ## For no sample size cap: if one arm is dropped then other arms 
      ## may receive more participants up to the nMax[sum(n_sample)] participants 
      if(capType=="noCap"){
        #print(pr.decision[,j-1])
        #drop_arm <- which(pr.decision[,j-1] < D0) + 1 # need to check
        drop_arm <- c(which(pr.decision[,j-1] < D0) + 1,which(pr.decision[,j-1] > D1) + 1)
        rand_pr[drop_arm] <- 0
        rand_pr <- rand_pr/sum(rand_pr)
        y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],gen_data_copd(N=N, bta=true_beta, sig=true_sig, rand_pr=rand_pr)[,3:4],SIMPLIFY=FALSE))
        n_store[j] <- nrow(y_init[[j]])
        n_interimSS[,j] <- table(y_init[[j]]$x)
        n_rand_pr[,j] <- rand_pr
      }
      else{
        y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],gen_data_copd(N=N, bta=true_beta, sig=true_sig, rand_pr=rand_pr)[,3:4],SIMPLIFY=FALSE))
        n_store[j] <- nrow(y_init[[j]])
        n_interimSS[,j] <- table(y_init[[j]]$x)
        n_rand_pr[,j] <- rand_pr
      }
      ##
      if(modelType=="Censored"){
        if(!is.null(L_per)){ # 
          #cen_n <- round(n_store[j]*L_per)
          cen_n <- round(n_sample[j]*L_per)
          LL <- sort(y_init[[j]]$y)[cen_n]
        }
        else{
          LL <- L
        }
        y_init[[j]]$y[y_init[[j]]$y<=LL] <- LL
        #per_missing[j] <- mean(y_init[[j]]$y<=LL)
        per_missing[,j] <- table(y_init[[j]][y_init[[j]]$y<=LL,"x"])/n_interimSS[,j] # armwise missing 
      }
      else{
        LL = NULL
        per_missing[,j] <- 0
      }
      ##
      if(unique(per_missing[,j])[length(unique(per_missing[,j]))]%in%0){
        model = model_noncen
        modelType2 = "nonCensored"
        LL = NULL
      }
      else{
        model = model_cen
        modelType2 = "Censored"
      }
      ##
      #print(c(j,modelType2,per_missing[j]))
      out[[j]] = trial_fast(data = y_init[[j]], 
                            model = model, # stan model from stan_model code 
                            modelType = modelType2, L=LL,
                            mu_beta = mu_beta, # with intercept -
                            sig_beta = sig_beta, # with intercept - 
                            a=a, b=b)
      ##
      pr.decision[,j] <- pnorm(0, out[[j]]$post.beta.mn[2:K], out[[j]]$post.beta.sd[2:K])
      ##
      ##
      ## for bino 
      y_init_bino[[j]] = y_init[[j]]
      y_init_bino[[j]][y_init[[j]]$y<=LL,"y"] = 1
      y_init_bino[[j]][y_init[[j]]$y>LL,"y"] = 0
      out_bino[[j]] = trial_fast_bernoulli(data = y_init_bino[[j]], 
                                           model = model_bino, # stan model from stan_model code 
                                           mu_beta = mu_beta, # with intercept -
                                           sig_beta = sig_beta)
      ck <- NULL
      for(l in 2:K){
        #ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))<1))
        ck <- c(ck,mean(exp(rnorm(10000,out_bino[[j]]$post.beta.mn[l],out_bino[[j]]$post.beta.sd[l]))>1))
      }
      pr.decision_bino[,j] <- ck
      ##
    }
    else{
      N <- n_sample[j]
      ## For no sample size cap: if one arm is dropped then other arms 
      ## may receive more participants up to the nMax[sum(n_sample)] participants 
      if(capType=="noCap"){
        #print(pr.decision[,j-1])
        #drop_arm <- which(pr.decision[,j-1] < D0) + 1
        drop_arm <- c(which(pr.decision[,j-1] < D0) + 1,which(pr.decision[,j-1] > D1) + 1)
        rand_pr[drop_arm] <- 0
        rand_pr <- rand_pr/sum(rand_pr)
        y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],gen_data_copd(N=N, bta=true_beta, sig=true_sig, rand_pr=rand_pr)[,3:4],SIMPLIFY=FALSE))
        n_store[j] <- nrow(y_init[[j]])
        n_interimSS[,j] <- table(y_init[[j]]$x)
        n_rand_pr[,j] <- rand_pr        
      }
      else{
        y_init[[j]] = data.frame(mapply(c,y_init[[j-1]],gen_data_copd(N=N, bta=true_beta, sig=true_sig, rand_pr=rand_pr)[,3:4],SIMPLIFY=FALSE))
        n_store[j] <- nrow(y_init[[j]])
        n_interimSS[,j] <- table(y_init[[j]]$x)
        n_rand_pr[,j] <- rand_pr
      }
      ##
      if(modelType=="Censored"){
        if(!is.null(L_per)){ # 
          #cen_n <- round((n_store[j]-n_store[j-1])*L_per)
          cen_n <- round(n_sample[j]*L_per)
          LL <- sort(y_init[[j]]$y)[cen_n]
        }
        else{
          LL <- L
        }
        y_init[[j]]$y[y_init[[j]]$y<=LL] <- LL
        #per_missing[j] <- mean(y_init[[j]]$y<=LL)
        per_missing[,j] <- table(y_init[[j]][y_init[[j]]$y<=LL,"x"])/n_interimSS[,j] # armwise missing 
      }
      else{
        LL = NULL
        per_missing[,j] <- 0
      }
      ##
      if(unique(per_missing[,j])[length(unique(per_missing[,j]))]%in%0){
        model = model_noncen
        modelType2 = "nonCensored"
        LL = NULL
      }
      else{
        model = model_cen
        modelType2 = "Censored"
      }
      ##
      #print(c(j,modelType2,per_missing[j]))
      out[[j]] = trial_fast(data = y_init[[j]], 
                            model = model, # stan model from stan_model code 
                            modelType = modelType2, L=LL,
                            mu_beta = mu_beta, # with intercept -
                            sig_beta = sig_beta, # with intercept - 
                            a=a, b=b)
      ##
      pr.decision[,j] <- pnorm(0, out[[j]]$post.beta.mn[2:K], out[[j]]$post.beta.sd[2:K])
      #drop_arm <- c(which(pr.decision[,j] < D0) + 1,which(pr.decision[,j] > D1) + 1)
      ##
      ##
      ## for bino 
      y_init_bino[[j]] = y_init[[j]]
      y_init_bino[[j]][y_init[[j]]$y<=LL,"y"] = 1
      y_init_bino[[j]][y_init[[j]]$y>LL,"y"] = 0
      out_bino[[j]] = trial_fast_bernoulli(data = y_init_bino[[j]], 
                                           model = model_bino, # stan model from stan_model code 
                                           mu_beta = mu_beta, # with intercept -
                                           sig_beta = sig_beta)
      ck <- NULL
      for(l in 2:K){
        #ck <- c(ck,mean(exp(rnorm(10000,out[[j]]$post.beta.mn[l],out[[j]]$post.beta.sd[l]))<1))
        ck <- c(ck,mean(exp(rnorm(10000,out_bino[[j]]$post.beta.mn[l],out_bino[[j]]$post.beta.sd[l]))>1))
      }
      pr.decision_bino[,j] <- ck
      ##
    }
  }
  ##
  post.beta.mn <- sapply(out, function(x) x$post.beta.mn)
  post.beta.sd <- sapply(out, function(x) x$post.beta.sd)
  post.sigma.mn <- sapply(out, function(x) x$post.sigma.mn)
  post.sigma.sd <- sapply(out, function(x) x$post.sigma.sd)
  ##
  dimnames(post.beta.mn)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(post.beta.sd)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(pr.decision)[[1]] = c(paste0("treat",1:(length(true_beta)-1)))
  ##
  dimnames(pr.decision_bino)[[1]] = c(paste0("treat",1:(length(true_beta)-1)))
  post.beta.mn_bino <- sapply(out_bino, function(x) x$post.beta.mn)
  post.beta.sd_bino <- sapply(out_bino, function(x) x$post.beta.sd)
  dimnames(post.beta.mn_bino)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(post.beta.sd_bino)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  ##
  dimnames(n_interimSS)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_interimSS)[[2]] = c(paste0("interim",1:n_j))
  ##
  dimnames(n_rand_pr)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(n_rand_pr)[[2]] = c(paste0("interim",1:n_j))
  ##
  per_missing <- per_missing*100
  dimnames(per_missing)[[1]] = c("control",paste0("treat",1:(length(true_beta)-1)))
  dimnames(per_missing)[[2]] = c(paste0("interim",1:n_j))
  #
  return(list(pr.decision=pr.decision,pr.decision_bino=pr.decision_bino,
              post.beta.mn=post.beta.mn, post.beta.sd=post.beta.sd,
              post.sigma.mn=post.sigma.mn, post.sigma.sd=post.sigma.sd,
              post.beta.mn_bino=post.beta.mn_bino, post.beta.sd_bino=post.beta.sd_bino,
              sample_size=n_interimSS, rand_pr=n_rand_pr, 
              n_store=n_store, per_missing=per_missing))
}
##
#true_beta <- c(1,-1,-0.5,-0.5)
#true_sig <- 1
#aa=trail_interim_copd_normal_bino(model_noncen=model_noncen, model_cen=model_cen, model_bino=model_bino,
#                               modelType = "Censored", L=log(0.016,base=10),L_per=NULL,     
#                               n_sample=c(100,100,100),
#                               rand_assignment=c(1,1,1,1),
#                               true_beta = true_beta,
#                               true_sig = true_sig,
#                               mu_beta = c(0,0,0,0), # with intercept -
#                               sig_beta = c(5,5,5,5), # with intercept - 
#                               a=3, b=1,
#                               capType="noCap", D0=0.1, D1=0.975)

##
## trial simulator
##
trial_interim_simulator_copd_normal_bino <- function(nSim=10, model_noncen=NULL, model_cen=NULL, model_bino=NULL,
                                    modelType = "Censored", L=log(0.016,base=10),L_per=NULL,     
                                    n_sample=c(100,50,50), # n_sample => increment values
                                    rand_assignment=c(1,1,1,1),
                                    true_beta = c(0.5,-1,-1,-1), # first input is intercept
                                    true_sig = c(1),
                                    prior_mu_beta = c(0,0,0,0), # with intercept -
                                    prior_sig_beta = sqrt(c(5,5,5,5)), # with intercept - 
                                    prior_IG_a=3, prior_IG_b=1,
                                    capType="noCap",D0=0.10,D1=0.975){
  ##
  #options(error = browser) 
  #options(error = recover)
  #options(error = NULL)
  library(pbapply)
  #out <- pblapply(1:nSim, 
  #                function(x) trail_interim_copd_normal_bino(model_noncen=model_noncen, model_cen=model_cen, model_bino=model_bino,
  #                                                                    modelType = modelType, L=L, L_per=L_per,
  #                                                                    n_sample=n_sample,
  #                                                                    rand_assignment=rand_assignment,
  #                                                                    true_beta = true_beta,
  #                                                                    true_sig = true_sig,
  #                                                                    mu_beta = prior_mu_beta, # with intercept -
  #                                                                    sig_beta = prior_sig_beta, # with intercept - 
  #                                                                    a=prior_IG_a, b=prior_IG_b,
  #                                                                    capType=capType, D0=D0, D1=D1))
  ##
  out <- pblapply(1:nSim, 
                  function(x) tryCatch(trail_interim_copd_normal_bino(model_noncen=model_noncen, model_cen=model_cen, model_bino=model_bino,
                                            modelType = modelType, L=L, L_per=L_per,
                                            n_sample=n_sample,
                                            rand_assignment=rand_assignment,
                                            true_beta = true_beta,
                                            true_sig = true_sig,
                                            mu_beta = prior_mu_beta, # with intercept -
                                            sig_beta = prior_sig_beta, # with intercept - 
                                            a=prior_IG_a, b=prior_IG_b,
                                            capType=capType, D0=D0, D1=D1), 
                                       error = function(e) NULL))
  
  ##
  out <- out[!sapply(out,is.null)]
  nSim <- length(out)
  ## keep at least 500 nSim to get the results
  if(nSim>=2500){
    ##
    #n_store <- out[[1]]$n_store/length(true_beta)
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
    tmp <- sapply(out, function(x) x$pr.decision_bino)
    pr.decision_bino <- array(tmp,dim=c(length(prior_mu_beta)-1,length(n_sample),nSim))
    dimnames(pr.decision_bino)[[1]] = paste0("treat",1:(length(prior_mu_beta)-1))
    dimnames(pr.decision_bino)[[2]] = paste0("interim",1:length(n_sample))
    ##
    tmp <- sapply(out, function(x) x$post.beta.mn_bino)
    post.beta.mn_bino <- array(tmp,dim=c(length(prior_mu_beta),length(n_sample),nSim))
    dimnames(post.beta.mn_bino)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
    dimnames(post.beta.mn_bino)[[2]] = paste0("interim",1:length(n_sample))
    ##
    tmp <- sapply(out, function(x) x$post.beta.sd_bino)
    post.beta.sd_bino <- array(tmp,dim=c(length(prior_mu_beta),length(n_sample),nSim))
    dimnames(post.beta.sd_bino)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
    dimnames(post.beta.sd_bino)[[2]] = paste0("interim",1:length(n_sample))
    ##
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
    post.sigma.mn <- sapply(out, function(x) x$post.sigma.mn)
    if(is.matrix(post.sigma.mn)){
      dimnames(post.sigma.mn)[[1]] = paste0("interim",1:length(n_sample))
    }
    else{
      post.sigma.mn <- as.matrix(post.sigma.mn)
      post.sigma.mn <- t(post.sigma.mn)
      dimnames(post.sigma.mn)[[1]] = paste0("interim",1:length(n_sample))
    }
    post.sigma.sd <- sapply(out, function(x) x$post.sigma.sd)
    if(is.matrix(post.sigma.sd)){
      dimnames(post.sigma.sd)[[1]] = paste0("interim",1:length(n_sample))
    }
    else{
      post.sigma.sd <- as.matrix(post.sigma.sd)
      post.sigma.sd <- t(post.sigma.sd)
      dimnames(post.sigma.sd)[[1]] = paste0("interim",1:length(n_sample))
    }
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
    per_missing <- sapply(out, function(x) x$per_missing)
    per_missing <- array(per_missing,dim=c(length(prior_mu_beta),length(n_sample),nSim))
    dimnames(per_missing)[[1]] = c("control",paste0("treat",1:(length(prior_mu_beta)-1)))
    dimnames(per_missing)[[2]] = paste0("interim",1:length(n_sample))
    #per_missing <- sapply(out, function(x) x$per_missing)
    #if(is.matrix(per_missing)){
    #  dimnames(per_missing)[[1]] = paste0("interim",1:length(n_sample))
    #}
    #else{
    #  per_missing <- as.matrix(per_missing)
    #  per_missing <- t(per_missing)
    #  dimnames(per_missing)[[1]] = paste0("interim",1:length(n_sample))
    #}
    ##
    return(list(sample_size=sample_size, 
                rand_pr=rand_pr,
                n_store=t(n_store),
                #per_missing=t(per_missing),
                per_missing=per_missing,
                pr.decision=pr.decision,
                pr.decision_bino=pr.decision_bino,
                post.beta.mn=post.beta.mn, 
                post.beta.sd=post.beta.sd,
                post.sigma.mn=t(post.sigma.mn), 
                post.sigma.sd=t(post.sigma.sd), 
                post.beta.mn_bino=post.beta.mn_bino, 
                post.beta.sd_bino=post.beta.sd_bino,
                n_sample=n_sample, nSim=nSim,
                true_beta=true_beta, true_sig=true_sig,
                prior_mu_beta=prior_mu_beta, prior_sig_beta=prior_sig_beta,
                prior_IG_a=prior_IG_a, prior_IG_b=prior_IG_b,
                capType=capType, modelType=modelType, 
                L=L, L_per=L_per, D0=D0, D1=D1))
  ##  
  }
  ##
  else{
    NULL
  }
}
##
##
fnc_replace_to_one <- function(x){
  ## x is a vector 
  if(length(which(x==1))>0){
    x[which(x==1)[1]:length(x)] <- 1
  }
  x
}

##

model_cen <- stan_model(file="lm_beta_leftcensored_optiz_knownL.stan")
model_noncen <- stan_model(file="lm_beta.stan")
model_bino <- stan_model(file="ber_logit.stan")

## simulation (#1)
n_sample=c(100,50,50)
rand_assignment=c(1,1,1,1)
true_bta = rbind(c(0.5,0,0,0),c(0.5,-0.5,-1,-1),c(0.5,0,-2,-1),
                 c(0,-0.5,-0.5,-1),c(0,-0.5,-1,-1),c(0,0,-2,-1),
                 c(-1,-0.5,-0.5,-1),c(-1,-0.5,-1,-1),c(-1,0,-2,-1))
true_sigma <- c(0.5,1,1.5)
D0 = 0.1; D1=0.975
nSim = 5000
id_scenario <- 0
store_res <- NULL
for(sig in 1:length(true_sigma)){
  for(bt in 1:nrow(true_bta)){
    res <- trial_interim_simulator_copd_normal_bino(nSim=nSim, model_noncen=model_noncen, 
                                             model_cen=model_cen, model_bino=model_bino,
                                             modelType = "Censored", L=log(0.016,base=10),L_per=NULL,     
                                             n_sample=n_sample, # n_sample => increment values
                                             rand_assignment=rand_assignment,
                                             true_beta = c(true_bta[bt,]), # first input is intercept
                                             true_sig = true_sigma[sig],
                                             prior_mu_beta = c(0,0,0,0), # with intercept -
                                             prior_sig_beta = sqrt(c(5,5,5,5)), # with intercept - 
                                             prior_IG_a=3, prior_IG_b=1,
                                             capType="noCap",D0=D0,D1=D1)
    ##
    ## summary stat
    ##
    dat <- data.table::CJ(
      interim = 1:length(n_sample),
      arm = 1:length(rand_assignment),
      id_scenario = NA_integer_,
      per_missing = NA_real_,
      per_missing_sd = NA_real_,
      true_beta = NA_real_,
      est_beta = NA_real_,
      est_beta_sd = NA_real_,
      true_sig = NA_real_,
      est_sig = NA_real_,
      est_sig_sd = NA_real_,
      est_bino_exp_beta = NA_real_,
      est_bino_exp_beta_sd = NA_real_,
      est_sample =  NA_real_,
      est_sd =  NA_real_
    )
    id_scenario <- id_scenario + 1
    dat$id_scenario <-  id_scenario
    dat$true_beta = rep(c(true_bta[bt,]),length(n_sample))
    dat$true_sig = true_sigma[sig]
    ##
    ## if no results are obtained from the OCs then we can't report any results 
    if(!is.null(res)){
      for(j in 1:length(n_sample)){
        #
        dat[dat$interim==j,"est_beta"] <- rowMeans(res$post.beta.mn[,j,])
        dat[dat$interim==j,"est_beta_sd"] <- rowMeans(res$post.beta.sd[,j,])
        #
        dat[dat$interim==j,"est_sig"] <-  mean(res$post.sigma.mn[,j])
        dat[dat$interim==j,"est_sig_sd"] <- mean(res$post.sigma.sd[,j])
        #
        dat[dat$interim==j,"est_bino_exp_beta"] <- rowMeans(exp(res$post.beta.mn_bino[,j,]))
        dat[dat$interim==j,"est_bino_exp_beta_sd"] <- apply(exp(res$post.beta.mn_bino[,j,]),1,sd)
        #
        dat[dat$interim==j,"est_sample"] <- rowMeans(res$sample_size[,j,])
        dat[dat$interim==j,"est_sd"] <- apply(res$sample_size[,j,],1,sd)
        #
        dat[dat$interim==j,"per_missing"] <- rowMeans(res$per_missing[,j,])
        dat[dat$interim==j,"per_missing_sd"] <- apply(res$per_missing[,j,],1,sd)
      }
      ##
      pr_decisionD0 <- NULL
      pr_decisionD1 <- NULL
      for(i in 1:(length(rand_assignment)-1)){
        tmp <- (res$pr.decision[i,,]<D0)
        pr_decisionD0 <- rbind(pr_decisionD0,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
        tmp <- (res$pr.decision[i,,]>D1)
        pr_decisionD1 <- rbind(pr_decisionD1,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
      }
      pr_decision <- cbind(rep(1:length(n_sample),(length(rand_assignment)-1)),c(t(pr_decisionD0)),c(t(pr_decisionD1)))
      dimnames(pr_decision)[[2]] <- c("interim","D0_0.1_Cont","D1_0.975_Cont")
      pr_decision <- data.frame(pr_decision)
      pr_decision$arm <- rep(2:length(rand_assignment),each=length(n_sample))
      ##
      dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
      ##
      pr_decisionD0 <- NULL
      pr_decisionD1 <- NULL
      for(i in 1:(length(rand_assignment)-1)){
        tmp <- (res$pr.decision_bino[i,,]<D0)
        pr_decisionD0 <- rbind(pr_decisionD0,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
        tmp <- (res$pr.decision_bino[i,,]>D1)
        pr_decisionD1 <- rbind(pr_decisionD1,rbind(rowMeans(apply(tmp,2,fnc_replace_to_one))))
      }
      pr_decision <- cbind(rep(1:length(n_sample),(length(rand_assignment)-1)),c(t(pr_decisionD0)),c(t(pr_decisionD1)))
      dimnames(pr_decision)[[2]] <- c("interim","D0_0.1_Bino","D1_0.975_Bino")
      pr_decision <- data.frame(pr_decision)
      pr_decision$arm <- rep(2:length(rand_assignment),each=length(n_sample))
      ##
      dat <- merge(dat,pr_decision,by=c("interim","arm"),all.x=TRUE)
      ##
    }
    ##
    else{
      dat$D0_0.1_Cont <- NA; dat$D0_0.1_Cont <- as.numeric(dat$D0_0.1_Cont)
      dat$D0_0.975_Cont <- NA; dat$D0_0.975_Cont <- as.numeric(dat$D0_0.975_Cont)
      dat$D0_0.1_Bino <- NA; dat$D0_0.1_Bino <- as.numeric(dat$D0_0.1_Bino)
      dat$D0_0.975_Bino <- NA; dat$D0_0.975_Bino <- as.numeric(dat$D0_0.975_Bino)
    }
    ##
    store_res <- rbind(store_res,dat, use.names=FALSE)
    write.csv(store_res,file=paste0("result_ContBino_D0.10_D1.975_",Sys.Date(),".csv"),row.names=FALSE)
    ##
  }
}

##

################################################################################
################################################################################
