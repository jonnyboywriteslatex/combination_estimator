Load required R libraries
```
library(boot)
library(quadprog)
```

## get\_est
Basic function to combine individual treatment effect estimates from trials and observational studies given in Section 2.2 of manuscript (B=number bootstrap iterations - if B=NA only the estimate will be returned)

Parameters:
* est\_rt: treatment effect estimate from randomized trial
* est\_os: treatment effect estimate from observational study
* sd\_rt: standard error of treatment effect estimate from randomized trial
* st\_os: standard error of treatment effect estimate from observational study
* B = number bootstrap resamples
* conf=confidence level
* logdata=TRUE (if the est_rt and est_os are given on log-scale (as might be true for logistic  regression)).  If logdata=TRUE, effects are exponentiated before outputing) 
``` 
get_est <- function(est_rt,est_os,sd_rt,sd_os,B=NA,conf=.95,logdata=TRUE){
  var_os <- sd_os^2
  var_rt <- sd_rt^2
  bias <- est_os-est_rt
  w <- (bias^2+var_os)/(bias^2+var_os+var_rt)
  est <- w*est_rt+(1-w)*est_os
  if(is.na(B)) return(est)
  # parametric bootstrap
  data <- c(est_rt=est_rt,est_os=est_os,sd_rt=sd_rt,sd_os=sd_os)
  est <- function(data){
    est_rt <- data[1]
    est_os <- data[2]
    var_rt <- data[3]^2
    var_os <- data[4]^2
    bias <- est_os-est_rt
    w <- (bias^2+var_os)/(bias^2+var_os+var_rt)
    return(w*est_rt+(1-w)*est_os)
  }
  the_est <- est(data)
  random_data <- function(data,mle=the_est){
    est_par <- mle[1]
    est_bias <- data[2]-data[1]
    sd_rt <- data[3]
    sd_os <- data[4]
    rt_b <- rnorm(1,mean=est_par,sd=sd_rt)
    os_b <- rnorm(1,mean=est_par+est_bias,sd=sd_os)
    return(c(rt_b,os_b,sd_rt,sd_os))
  }
  the.boot <- boot(data, est, R = B, sim ="parametric", ran.gen = random_data, mle = mean(the_est))
  theints <- boot.ci(the.boot,conf=conf,type=c("norm","basic", "stud", "perc"))
  if(logdata) return(c(exp(the_est),exp(theints$norm[2:3])))
 return(c(the_est[1],theints$norm[2:3]))
}
```
##  est\_cons 
Function to calculate basic estimate given in Section 2.5  of manuscript tapering the weight for the trial toward 1.

Parameters:
* est\_rt: treatment effect estimate from randomized trial
* est\_os: treatment effect estimate from observational study
* sd\_rt: standard error of treatment effect estimate from randomized trial
* st\_os: standard error of treatment effect estimate from observational study
* B = number bootstrap resamples
* conf=confidence level* logdata=TRUE (if the est_rt and est_os are given on log-scale (as might be true for logistic  regression)).  If logdata=TRUE, effects are exponentiated before outputing) 
``` 
est_cons <- function(est_rt,est_os,sd_rt,sd_os,B=NA,conf=.95,logdata=TRUE){
  if(is.na(B)){
    var_os <- sd_os^2
    var_rt <- sd_rt^2
    bias <- est_os-est_rt
    bias <-max(abs(bias-2*sqrt(var_rt+var_os)),abs(bias+2*sqrt(var_rt+var_os)))
    w <- (bias^2+var_os)/(bias^2+var_os+var_rt)
    return(w*est_rt+(1-w)*est_os)
  }
  data <- c(est_rt=est_rt,est_os=est_os,sd_rt=sd_rt,sd_os=sd_os)
  est <- function(data){
    est_rt <- data[1]
    est_os <- data[2]
    var_rt <- data[3]^2
    var_os <- data[4]^2
    bias <- est_os-est_rt
    bias <- max(abs(bias-2*sqrt(var_rt+var_os)),abs(bias+2*sqrt(var_rt+var_os)))
    w <-  (bias^2+var_os)/(bias^2+var_os+var_rt)
    return(w*est_rt+(1-w)*est_os)
  }
  the_est <- est(data)
random_data <- function(data,mle=the_est){
est_par <- mle[1]
    est_bias <- data[2]-data[1]
    sd_rt <- data[3]
    sd_os <- data[4]
    rt_b <- rnorm(1,mean=est_par,sd=sd_rt)
    os_b <- rnorm(1,mean=est_par+est_bias,sd=sd_os)
    return(c(rt_b,os_b,sd_rt,sd_os))
  }
  the.boot <- boot(data, est, R = B, sim =
"parametric",
                   ran.gen = random_data, mle =mean(the_est))
  theints <- boot.ci(the.boot,conf=conf,type=c("norm","basic","stud", "perc"))
  if(logdata) return(c(exp(the_est),exp(theints$norm[2:3])))
  return(c(the_est[1],theints$norm[2:3]))
}
```
##  fe_meta
Fixed effects meta analysis

Parameters
* theta\_vec,se\_vec both vectors
 
```
fe_meta <- function(theta_vec,se_vec){
 weights <- 1/se_vec^2/(sum(1/se_vec^2))
 return(c(sum(weights*theta_vec),sum(weights^2*se_vec^2)))
}
```

##  re_meta
Random effects meta analysis using DerSimmion/Laird estimate for between study variance
 
Parameters
* theta\_vec,se\_vec both vectors

```
re_meta <- function(theta_vec,se_vec){
  F <- fe_meta(theta_vec,se_vec)[1]
  Q <- sum(se_vec^-2*(theta_vec-F)^2)
  sigmab2 <- max(0, (Q-(length(theta_vec)-1))/(sum(se_vec^-2)-sum(se_vec^-4)/sum(se_vec^-2)))
  weights <- 1/(se_vec^2+sigmab2)/(sum(1/(se_vec^2+sigmab2)))
return(c(sum(weights*theta_vec),sum(weights^2*(se_vec^2+sigmab2))))
}
```

##  est_fixedeffect
Fixed-bias meta analysis (from 3.1).  The estimate for the randomized trials is calculated using a fixed effects meta analysis. 
This function allows negative weights for the observational studies

Parameters:
* est\_rt: treatment effect vector from collection of randomized trials 
* est\_os: treatment effect vector from collection of observational studies 
* sd\_rt: standard error of treatment effect estimates from randomized trial
* st\_os: standard error of treatment effect estimates from observational study
* B = number bootstrap resamples
* logdata=TRUE (if the est_rt and est_os are given on log-scale (as might be true for logistic  regression).  If logdata=TRUE, effects are exponentiated before outputing)

```
est_fixedeffect <- function(est_rt,est_os,sd_rt,sd_os,B=1000,logdata=TRUE,conf=.95){
  if(is.na(B)){
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat
    var_os <- sd_os^2
    var_rt <- sd_rt^2
    Dmat <- 1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))
    bvec <- c(1)
    meq <- 1
    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))
    Amat <- t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))
    weights <- solve.QP(Dmat=Dmat,dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution
    w_rt <- weights[1:length(est_rt)]
    w_os <- weights[(length(est_rt)+1):length(weights)]
    if(logdata) return(exp((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))))
   return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))
  }
  N_rt <- length(est_rt)
  N_os <- length(est_os)
  data <- c(est_rt,sd_rt,est_os,sd_os,N_rt,N_os)
  get_est <- function(data){
    N_os <- data[length(data)]
    N_rt <- data[length(data)-1]
    est_rt <- data[1:N_rt]
    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]
    sd_rt <- data[(N_rt+1):(2*N_rt)]
    sd_os <- data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat
    var_os <- sd_os^2
    var_rt <- sd_rt^2
    Dmat <- 1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))
    bvec <- c(1)
    meq <- 1
    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))
    Amat <- t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))
    weights <- solve.QP(Dmat=Dmat,dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution
    w_rt <- weights[1:length(est_rt)]
    w_os <- weights[(length(est_rt)+1):length(weights)]
    return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))
  }
  the_est <- get_est(data)
random_data <- function(data,mle=the_est){
    N_os <- data[length(data)]
    N_rt <- data[length(data)-1]
    est_rt <- data[1:N_rt]
    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]
    sd_rt <- data[(N_rt+1):(2*N_rt)]
    sd_os <- data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat
    rt_b <- rnorm(length(sd_rt),mean=mle,sd=sd_rt)
    os_b <- rnorm(length(sd_os),mean=mle+est_bias,sd=sd_os)
    return(c(rt_b,sd_rt,os_b,sd_os,N_rt,N_os))
  }

  the.boot <- boot(data, get_est, R = B, sim = "parametric", ran.gen = random_data, mle = the_est)
  theints <- boot.ci(the.boot,conf=conf,type=c("norm","basic", "stud", "perc")) 
  if(logdata) return(c(exp(the_est),exp(theints$norm[2:3])))
  return(c(the_est,theints$norm[2:3]))
}
```

##  est_fixedeffect_nn
Fixed-bias meta analysis (from 3.1).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  This function does not allow negative weights for the observational studies


Parameters:
* est\_rt: treatment effect vector from collection of randomized trials 
* est\_os: treatment effect vector from collection of observational studies 
* sd\_rt: standard error of treatment effect estimates from randomized trial
* st\_os: standard error of treatment effect estimates from observational study
* B = number bootstrap resamples
* logdata=TRUE (if the est_rt and est_os are given on log-scale (as might be true for logistic  regression).  If logdata=TRUE, effects are exponentiated before outputing) 

``` 
est_fixedeffect_nn <- function(est_rt,est_os,sd_rt,sd_os,B=1000,logdata=TRUE,conf=.95){
  if(is.na(B)){
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat
    var_os <- sd_os^2
    var_rt <- sd_rt^2
    Dmat <- 1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))
    bvec <- c(1,rep(0,length(est_rt)),rep(0,length(est_os)))
    meq <- 1
    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))
    Amat <- t(matrix(rbind(c(rep(1,length(est_rt)),rep(1,length(est_os))),diag(1,length(est_rt)+length(est_os))),nrow=length(bvec)))
    weights <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution
    w_rt <- weights[1:length(est_rt)]
    w_os <- weights[(length(est_rt)+1):length(weights)]
       if(logdata) return(exp((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))))
   return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))
  }
  N_rt <- length(est_rt)
  N_os <- length(est_os)
  data <- c(est_rt,sd_rt,est_os,sd_os,N_rt,N_os)
  get_est <- function(data){
    N_os <- data[length(data)]
    N_rt <- data[length(data)-1]
    est_rt <- data[1:N_rt]
    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]
    sd_rt <- data[(N_rt+1):(2*N_rt)]
    sd_os <- data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat
    var_os <- sd_os^2
    var_rt <- sd_rt^2
    Dmat <- 1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))
    bvec <- c(1,rep(0,length(est_rt)),rep(0,length(est_os)))
    meq <- 1
    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))
    Amat <- t(matrix(rbind( c(rep(1,length(est_rt)),rep(1,length(est_os))),diag(1,length(est_rt)+length(est_os))),nrow=length(bvec)))
    weights <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution
    w_rt <- weights[1:length(est_rt)]
    w_os <- weights[(length(est_rt)+1):length(weights)]
return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))) 
  }
  the_est <- get_est(data)
  random_data <- function(data,mle=the_est){  
    N_os <- data[length(data)]
    N_rt <- data[length(data)-1]
    est_rt <- data[1:N_rt]
    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]
    sd_rt <- data[(N_rt+1):(2*N_rt)]
    sd_os <- data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat
    rt_b <- rnorm(length(sd_rt),mean=mle,sd=sd_rt)
    os_b <- rnorm(length(sd_os),mean=mle+est_bias,sd=sd_os)
    return(c(rt_b,sd_rt,os_b,sd_os,N_rt,N_os))
  }
  the.boot <- boot(data, get_est, R = B, sim = "parametric", ran.gen = random_data, mle = the_est)
  theints <- boot.ci(the.boot,conf=conf,type=c("norm","basic", "stud", "perc")) 
  if(logdata) return(c(exp(the_est),exp(theints$norm[2:3])))
  return(c(the_est,theints$norm[2:3]))
}
```

##  est_randomeffect
Random-bias meta analysis (from Section 3.2 of manuscript) using DerSimmion/Laird estimate for between study variance (theta\_vec,se\_vec are vectors).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  

Parameters
* est\_rt: treatment effect vector from collection of randomized trials
* est\_os: treatment effect vector from collection of observational studies
* sd\_rt: standard error of treatment effect estimates from randomized trial
* st\_os: standard error of treatment effect estimates from observational study
* B = number bootstrap resamples
* logdata=TRUE (if the est_rt and est_os are given on log-scale (as might be true for logistic  regression).  If logdata=TRUE, effects are exponentiated before outputing) 
```
est_randomeffect <- function(est_rt,est_os,sd_rt,sd_os,B=1000,logdata=TRUE,conf=.95){
  if(is.na(B)){
    N_os <- length(est_os)
    N_rt <- length(est_rt)
    out <- fe_meta(est_rt,sd_rt)
    theta_hat <- out[1]
    theta_hat_os <-  fe_meta(est_os,sd_os)[1]
    Q <- sum(sd_os^-2*(est_os-theta_hat_os)^2)
    sigmab2 <- max(0, (Q-(length(est_os)-1))/(sum(sd_os^-2)-sum(sd_os^-4)/sum(sd_os^-2)))
    est_bias <- rep(fe_meta(est_os,sqrt(sd_os^2+sigmab2))[1] - theta_hat,length(sd_os))
    var_os <- sd_os^2 + sigmab2
    var_rt <- sd_rt^2
    Dmat <- 1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,var_os)))
    bvec <- c(1)
    meq <- 1
    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))
    Amat <- t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))
    weights <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution
    w_rt <- weights[1:length(est_rt)]
    w_os <- weights[(length(est_rt)+1):length(weights)]
    if(logdata) return(exp((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))))
return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))
  }
    N_rt <- length(est_rt)
  N_os <- length(est_os)
  data <- c(est_rt,sd_rt,est_os,sd_os,N_rt,N_os)
  est <- function(data){
    N_os <- data[length(data)]
    N_rt <- data[length(data)-1]
    est_rt <- data[1:N_rt]
    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]
    sd_rt <- data[(N_rt+1):(2*N_rt)]
    sd_os <- data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]
    out <- fe_meta(est_rt,sd_rt)
    theta_hat <- out[1]
    theta_hat_os <-  fe_meta(est_os,sd_os)[1]
    Q <- sum(sd_os^-2*(est_os-theta_hat_os)^2)
    sigmab2 <- max(0, (Q-(length(est_os)-1))/(sum(sd_os^-2)-sum(sd_os^-4)/sum(sd_os^-2)))
    est_bias <- rep(fe_meta(est_os,sqrt(sd_os^2+sigmab2))[1] - theta_hat,length(sd_os))
    var_os <- sd_os^2 + sigmab2
    var_rt <- sd_rt^2
    Dmat <- 1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,var_os)))
    bvec <- c(1)
    meq <- 1
    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))
    Amat <- t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))
    weights <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution
    w_rt <- weights[1:length(est_rt)]
    w_os <- weights[(length(est_rt)+1):length(weights)]
   return(c((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)),sigmab2
  ))
}
the_est <- est(data)
  random_data <- function(data,mle=the_est){
    N_os <- data[length(data)]
    N_rt <- data[length(data)-1]
    est_rt <- data[1:N_rt]
    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]
    sd_rt <- data[(N_rt+1):(2*N_rt)]
    sd_os <- data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]
    theta_hat <- fe_meta(est_rt,sd_rt)[1]
    est_bias <- est_os-theta_hat+rnorm(N_os,mean=0,sd=sqrt(mle[2]))
    rt_b <- rnorm(length(sd_rt),mean=mle[1],sd=sd_rt)
    os_b <-rnorm(length(sd_os),mean=mle[1]+est_bias,sd=sd_os)
    return(c(rt_b,sd_rt,os_b,sd_os,N_rt,N_os))
  }
  the.boot <- boot(data, est, R = B, sim ="parametric",ran.gen = random_data, mle = the_est)
  theints <-boot.ci(the.boot,conf=conf,type=c("norm","basic","stud", "perc")) 
  if(logdata) return(c(exp(the_est[1]),exp(theints$norm[2:3])))
  return(c(the_est[1],theints$norm[2:3]))
}
```
##  Examples

```

###  some examples - from section 4.2

# function to get log-scale standard errors and estimates from OR confidence limits
getse <- function(l1,l2){
 log1 <- log(l1)
  log2 <- log(l2)
  a <- ((log1+log2)/2)
  sel <- (log2-log1)/(2*1.96)
  return(c(a,sel))
}
## estimates from randomized trial
rtmat <- matrix(0,ncol=2,nrow=2)
rtmat[1,] <- getse(0.25,1) 
rtmat[2,] <- getse(0.42,0.97) 

## estimates from observational study
osmat <- matrix(0,ncol=2,nrow=19)
osmat[1,] <- getse(0.3,0.64) 
osmat[2,] <- getse(0.35,0.86) 
osmat[3,] <- getse(0.54,0.78) 
osmat[4,] <- getse(0.18,1.01) 
osmat[5,] <- getse(0.55,0.85) 
osmat[6,] <- getse(0.49,1.09) 
osmat[7,] <- getse(0.42,0.79) 
osmat[8,] <- getse(0.40,1.15) 
osmat[9,] <- getse(0.3,0.87) 
osmat[10,] <- getse(0.4,0.71) 
osmat[11,] <- getse(0.25,1.91) 
osmat[12,] <- getse(0.61,0.95) 
osmat[13,] <- getse(0.51,1.23) 
osmat[14,] <- getse(0.34,1.14) 
osmat[15,] <- getse(0.20,.45) 
osmat[16,] <- getse(0.19,1.05) 
osmat[17,] <- getse(0.15,.59) 
osmat[18,] <- getse(0.19,0.79) 
osmat[19,] <- getse(0.12,0.52) 
 
# Section 3.1 approach
# not restricting weights
est_fixedeffect(rtmat[,1],osmat[,1],rtmat[,2],osmat[,2],B=10000)
# restricting weights to be negative
est_fixedeffect_nn(rtmat[,1],osmat[,1],rtmat[,2],osmat[,2],B=10000)

# random effects meta analysis
est_randomeffect(rtmat[,1],osmat[,1],rtmat[,2],osmat[,2],B=10000)

# Separate meta analyses and then combining using methods from Section 2.1 (since estimated sigma_b > 0 use random effects meta analysis for observational studies.
stuff1 <- fe_meta(rtmat[,1],rtmat[,2])
stuff2 <- re_meta(osmat[,1],osmat[,2])

get_est(est_rt=stuff1[1],est_os=stuff2[1],sd_rt=sqrt(stuff1[2]),sd_os=sqrt(stuff2[2]),B=10000)

# note that the estimate is equal to that produced using est_randomeffect
```
