
library(boot)
library(quadprog)

# Basic function to combine randomized trials and observational data (B=number bootstrap iterations - if B=NA only the estimate will be returned)

##  get\_est - function to calculate basic estimate given in Section 2.2  of manuscript.

##  parameters:

##  est\_rt: treatment effect estimate from randomized trial

##  est\_os: treatment effect estimate from observational study

##  sd\_rt: standard error of treatment effect estimate from randomized trial

##  st\_os: standard error of treatment effect estimate from observational study

##  B = number bootstrap resamples
 
```
get_est
<- function(est_rt,est_os,sd_rt,sd_os,B=NA,conf=.95){

  var_os <- sd_os^2

  var_rt <- sd_rt^2

  bias <- est_os-est_rt

  w <- 
(bias^2+var_os)/(bias^2+var_os+var_rt)

  est <- w*est_rt+(1-w)*est_os

  

  if(is.na(B)) return(est)

  # parametric bootstrap

  data <-
c(est_rt=est_rt,est_os=est_os,sd_rt=sd_rt,sd_os=sd_os)

  est <- function(data){

    est_rt <- data[1]

    est_os <- data[2]

    var_rt <- data[3]^2

    var_os <- data[4]^2

    bias <- est_os-est_rt

    w <- 
(bias^2+var_os)/(bias^2+var_os+var_rt)

    return(w*est_rt+(1-w)*est_os)

    

  }

  the_est <- est(data)

  

  random_data <- function(data,mle=the_est){

    est_par <- mle[1]

    est_bias <- data[2]-data[1]

    sd_rt <- data[3]

    sd_os <- data[4]

    rt_b <- rnorm(1,mean=est_par,sd=sd_rt)

    os_b <-
rnorm(1,mean=est_par+est_bias,sd=sd_os)

    return(c(rt_b,os_b,sd_rt,sd_os))

  }

  the.boot <- boot(data, est, R = B, sim =
"parametric",

                   ran.gen = random_data, mle =
mean(the_est))

  theints <-
boot.ci(the.boot,conf=conf,type=c("norm","basic",
"stud", "perc"))

  return(c(the_est[1],theints$norm[2:3]))

  

}
```
 

##  est\_cons - function to calculate basic estimate given in Section 2.2  of manuscript tapering the weight for the trial toward 1.
##  parameters:
##  est\_rt: treatment effect estimate from randomized trial
##  est\_os: treatment effect estimate from observational study
##  sd\_rt: standard error of treatment effect estimate from randomized trial
##  st\_os: standard error of treatment effect estimate from observational study
##  B = number bootstrap resamples

```
est_cons
<- function(est_rt,est_os,sd_rt,sd_os,B=NA,conf=.95){

  if(is.na(B)){

    var_os <- sd_os^2

    var_rt <- sd_rt^2

    bias <- est_os-est_rt

    bias <-
max(abs(bias-2*sqrt(var_rt+var_os)),abs(bias+2*sqrt(var_rt+var_os)))

    w <- 
(bias^2+var_os)/(bias^2+var_os+var_rt)

    return(w*est_rt+(1-w)*est_os)

  }

  data <-
c(est_rt=est_rt,est_os=est_os,sd_rt=sd_rt,sd_os=sd_os)

  est <- function(data){

    est_rt <- data[1]

    est_os <- data[2]

    var_rt <- data[3]^2

    var_os <- data[4]^2

    bias <- est_os-est_rt

    bias <-
max(abs(bias-2*sqrt(var_rt+var_os)),abs(bias+2*sqrt(var_rt+var_os)))

    w <- 
(bias^2+var_os)/(bias^2+var_os+var_rt)

    return(w*est_rt+(1-w)*est_os)

    

  }

  the_est <- est(data)

  

  random_data <- function(data,mle=the_est){

    est_par <- mle[1]

    est_bias <- data[2]-data[1]

    sd_rt <- data[3]

    sd_os <- data[4]

    rt_b <- rnorm(1,mean=est_par,sd=sd_rt)

    os_b <-
rnorm(1,mean=est_par+est_bias,sd=sd_os)

    return(c(rt_b,os_b,sd_rt,sd_os))

  }

  the.boot <- boot(data, est, R = B, sim =
"parametric",

                   ran.gen = random_data, mle =
mean(the_est))

  theints <-
boot.ci(the.boot,conf=conf,type=c("norm","basic",
"stud", "perc"))



  return(c(the_est[1],theints$norm[2:3]))

  

}

 
```
#  fixed effects meta analysis (theta\_vec,se\_vec are vectors)
 
```
fe_meta
<- function(theta_vec,se_vec){

  

  weights <- 1/se_vec^2/(sum(1/se_vec^2))

 
return(c(sum(weights*theta_vec),sum(weights^2*se_vec^2)))

  

}```

 #  random effects meta analysis using DerSimmion/Laird estimate for between study variance (theta\_vec,se\_vec are vectors)
 
```
re_meta
<- function(theta_vec,se_vec){

  F <- fe_meta(theta_vec,se_vec)[1]

  Q <- sum(se_vec^-2*(theta_vec-F)^2)

  sigmab2 <- max(0,
(Q-(length(theta_vec)-1))/(sum(se_vec^-2)-sum(se_vec^-4)/sum(se_vec^-2)))

  weights <-
1/(se_vec^2+sigmab2)/(sum(1/(se_vec^2+sigmab2)))

 
return(c(sum(weights*theta_vec),sum(weights^2*(se_vec^2+sigmab2))))

}

 

 ```

#  Fixed-bias meta analysis (from 3.1).  The estimate for the randomized trials is
calculated using a fixed effects meta analysis. 
This function allows negative weights for the observational studies

##  est_rt: treatment effect vector from
collection of randomized trials

##  est_os: treatment effect vector from
collection of observational studies

##  sd_rt: standard error of treatment effects
from randomized trials

##  st_os: standard error of treatment effects
from observational studies

##  B = number bootstrap resamples

##  logdata=TRUE (if the est_rt and est_os are
given on log-scale (as might be true for logistic regression).  If logdata=TRUE, effects are exponentiated
before outputing

 
```
est_fixedeffect
<- function(est_rt,est_os,sd_rt,sd_os,B=1000,logdata=TRUE,conf=.95){

  if(is.na(B)){

    theta_hat <- fe_meta(est_rt,sd_rt)[1]

    est_bias <- est_os-theta_hat

    var_os <- sd_os^2

    var_rt <- sd_rt^2

    Dmat <-
1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))

    bvec <- c(1)

    meq <- 1

    dvec <-
c(rep(0,length(est_rt)),rep(0,length(est_os)))

    Amat <-
t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))

    weights <- solve.QP(Dmat=Dmat,
dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w_rt <- weights[1:length(est_rt)]

    w_os <-
weights[(length(est_rt)+1):length(weights)]

    if(logdata)
return(exp((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))))

   
return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))

  }

  N_rt
<- length(est_rt)

  N_os <- length(est_os)

  data <-
c(est_rt,sd_rt,est_os,sd_os,N_rt,N_os)

  get_est <- function(data){

    N_os <- data[length(data)]

    N_rt <- data[length(data)-1]

    est_rt <- data[1:N_rt]

    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]

    sd_rt <- data[(N_rt+1):(2*N_rt)]

    sd_os <-
data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]

    theta_hat <- fe_meta(est_rt,sd_rt)[1]

    est_bias <- est_os-theta_hat

    var_os <- sd_os^2

    var_rt <- sd_rt^2

    Dmat <-
1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))

    bvec <- c(1)

    meq <- 1

    dvec <-
c(rep(0,length(est_rt)),rep(0,length(est_os)))

    Amat <-
t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))

    weights <- solve.QP(Dmat=Dmat,
dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w_rt <- weights[1:length(est_rt)]

    w_os <-
weights[(length(est_rt)+1):length(weights)]

    return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))


  }

  

  the_est <- get_est(data)

  

  random_data <- function(data,mle=the_est){

    

    

    N_os <- data[length(data)]

    N_rt <- data[length(data)-1]

    est_rt <- data[1:N_rt]

    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]

    sd_rt <- data[(N_rt+1):(2*N_rt)]

    sd_os <-
data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]

    theta_hat <- fe_meta(est_rt,sd_rt)[1]

    est_bias <- est_os-theta_hat

    rt_b <-
rnorm(length(sd_rt),mean=mle,sd=sd_rt)

    os_b <- rnorm(length(sd_os),mean=mle+est_bias,sd=sd_os)

    return(c(rt_b,sd_rt,os_b,sd_os,N_rt,N_os))

  }

  the.boot <- boot(data, get_est, R = B, sim
= "parametric",

                   ran.gen = random_data, mle =
the_est)

  theints <-
boot.ci(the.boot,conf=conf,type=c("norm","basic",
"stud", "perc")) 

  if(logdata)
return(c(exp(the_est),exp(theints$norm[2:3])))

  return(c(the_est,theints$norm[2:3]))

}
```

#  Fixed-bias meta analysis (from 3.1).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  This function does not allow negative weights for the observational studies
##  est\_rt: treatment effect vector from collection of randomized trials
##  est\_os: treatment effect vector from collection of observational studies
##  sd\_rt: standard error of treatment effects from randomized trials
##  st\_os: standard error of treatment effects from observational studies
##  B = number bootstrap resamples
##  logdata=TRUE (if the est\_rt and est\_os are given on log-scale (as might be true for logistic regression).  If logdata=TRUE, effects are exponentiated before outputing


```
est_fixedeffect_nn
<- function(est_rt,est_os,sd_rt,sd_os,B=1000,logdata=TRUE,conf=.95){

  if(is.na(B)){

    theta_hat <- fe_meta(est_rt,sd_rt)[1]

    est_bias <- est_os-theta_hat

    var_os <- sd_os^2

    var_rt <- sd_rt^2

    Dmat <-
1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))

    bvec <-
c(1,rep(0,length(est_rt)),rep(0,length(est_os)))

    meq <- 1

    dvec <-
c(rep(0,length(est_rt)),rep(0,length(est_os)))

    Amat <- t(matrix(rbind(
c(rep(1,length(est_rt)),rep(1,length(est_os))),diag(1,length(est_rt)+length(est_os))),nrow=length(bvec)))

    weights <- solve.QP(Dmat=Dmat, dvec=dvec,
Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w_rt <- weights[1:length(est_rt)]

    w_os <-
weights[(length(est_rt)+1):length(weights)]

    

    if(logdata)
return(exp((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))))

   
return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))

  }

  N_rt <- length(est_rt)

  N_os <- length(est_os)

  data <-
c(est_rt,sd_rt,est_os,sd_os,N_rt,N_os)

  get_est <- function(data){

    N_os <- data[length(data)]

    N_rt <- data[length(data)-1]

    est_rt <- data[1:N_rt]

    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]

    sd_rt <- data[(N_rt+1):(2*N_rt)]

    sd_os <-
data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]

    theta_hat <- fe_meta(est_rt,sd_rt)[1]

    est_bias <- est_os-theta_hat

    var_os <- sd_os^2

    var_rt <- sd_rt^2

    Dmat <-
1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,sd_os^2)))

    bvec <-
c(1,rep(0,length(est_rt)),rep(0,length(est_os)))

    meq <- 1

    dvec <- c(rep(0,length(est_rt)),rep(0,length(est_os)))

    Amat <- t(matrix(rbind(
c(rep(1,length(est_rt)),rep(1,length(est_os))),diag(1,length(est_rt)+length(est_os))),nrow=length(bvec)))

    weights <- solve.QP(Dmat=Dmat,
dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w_rt <- weights[1:length(est_rt)]

    w_os <-
weights[(length(est_rt)+1):length(weights)]

   
return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))) 

  }

  

  the_est <- get_est(data)

  

  random_data <- function(data,mle=the_est){

    

    

    N_os <- data[length(data)]

    N_rt <- data[length(data)-1]

    est_rt <- data[1:N_rt]

    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]

    sd_rt <- data[(N_rt+1):(2*N_rt)]

    sd_os <-
data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]

    theta_hat <- fe_meta(est_rt,sd_rt)[1]

    est_bias <- est_os-theta_hat

    rt_b <-
rnorm(length(sd_rt),mean=mle,sd=sd_rt)

    os_b <-
rnorm(length(sd_os),mean=mle+est_bias,sd=sd_os)

    return(c(rt_b,sd_rt,os_b,sd_os,N_rt,N_os))

  }

  the.boot <- boot(data, get_est, R = B, sim
= "parametric",

                   ran.gen = random_data, mle =
the_est)

  theints <-
boot.ci(the.boot,conf=conf,type=c("norm","basic",
"stud", "perc")) 

  if(logdata)
return(c(exp(the_est),exp(theints$norm[2:3])))

  return(c(the_est,theints$norm[2:3]))

}

 ```

#  Random-bias meta analysis (from 3.1) using DerSimmion/Laird estimate for between study variance (theta\_vec,se\_vec are vectors).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  This function does not allow negative weights for the trials
##  est\_rt: treatment effect vector from collection of randomized trials
##  est\_os: treatment effect vector from collection of observational studies
##  sd\_rt: standard error of treatment effects from randomized trials
##  st\_os: standard error of treatment effects from observational studies
##  B = number bootstrap resamples
##  logdata=TRUE (if the est\_rt and est\_os are given on log-scale (as might be true for logistic regression).  If logdata=TRUE, effects are exponentiated before outputing

```
est_randomeffect
<- function(est_rt,est_os,sd_rt,sd_os,B=1000,logdata=TRUE,conf=.95){

  

  if(is.na(B)){

    N_os <- length(est_os)

    N_rt <- length(est_rt)

    out <- fe_meta(est_rt,sd_rt)

    theta_hat <- out[1]

    theta_hat_os <-  fe_meta(est_os,sd_os)[1]

    Q <-
sum(sd_os^-2*(est_os-theta_hat_os)^2)

    sigmab2 <- max(0,
(Q-(length(est_os)-1))/(sum(sd_os^-2)-sum(sd_os^-4)/sum(sd_os^-2)))

    est_bias <-
rep(fe_meta(est_os,sqrt(sd_os^2+sigmab2))[1] - theta_hat,length(sd_os))

    var_os <- sd_os^2 + sigmab2

    var_rt <- sd_rt^2

    Dmat <-
1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,var_os)))

    bvec <- c(1)

    meq <- 1

    dvec <-
c(rep(0,length(est_rt)),rep(0,length(est_os)))

    Amat <-
t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))

    weights <- solve.QP(Dmat=Dmat,
dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w_rt <- weights[1:length(est_rt)]

    w_os <- weights[(length(est_rt)+1):length(weights)]

    if(logdata)
return(exp((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt))))

   
return((sum(w_rt*est_rt)+sum(w_os*est_os))/(sum(w_os)+sum(w_rt)))

  }

  

  N_rt <- length(est_rt)

  N_os <- length(est_os)

  data <-
c(est_rt,sd_rt,est_os,sd_os,N_rt,N_os)

  est <- function(data){

    N_os <- data[length(data)]

    N_rt <- data[length(data)-1]

    est_rt <- data[1:N_rt]

    est_os <- data[(2*N_rt+1):(2*N_rt+N_os)]

    sd_rt <- data[(N_rt+1):(2*N_rt)]

    sd_os <-
data[(2*N_rt+N_os+1):(2*N_rt+2*N_os)]

    out <- fe_meta(est_rt,sd_rt)

    theta_hat <- out[1]

    theta_hat_os <-  fe_meta(est_os,sd_os)[1]

    Q <-
sum(sd_os^-2*(est_os-theta_hat_os)^2)

    sigmab2 <- max(0,
(Q-(length(est_os)-1))/(sum(sd_os^-2)-sum(sd_os^-4)/sum(sd_os^-2)))

    est_bias <-
rep(fe_meta(est_os,sqrt(sd_os^2+sigmab2))[1] - theta_hat,length(sd_os))

    var_os <- sd_os^2 + sigmab2

    var_rt <- sd_rt^2

    Dmat <-
1*(c(rep(0,length(est_rt)),est_bias)%*%t(c(rep(0,length(est_rt)),est_bias))+diag(c(sd_rt^2,var_os)))

    bvec <- c(1)

    meq <- 1

    dvec <-
c(rep(0,length(est_rt)),rep(0,length(est_os)))

    Amat <-
t(matrix(c(rep(1,length(est_rt)),rep(1,length(est_os))),nrow=length(bvec)))

    weights <- solve.QP(Dmat=Dmat,
dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w_rt <- weights[1:length(est_rt)]

    w_os <-
weights[(length(est_rt)+1):length(weights)]

   
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

    est_bias <-
est_os-theta_hat+rnorm(N_os,mean=0,sd=sqrt(mle[2]))

    rt_b <-
rnorm(length(sd_rt),mean=mle[1],sd=sd_rt)

    os_b <-
rnorm(length(sd_os),mean=mle[1]+est_bias,sd=sd_os)

    return(c(rt_b,sd_rt,os_b,sd_os,N_rt,N_os))

  }

  the.boot <- boot(data, est, R = B, sim =
"parametric",

                   ran.gen = random_data, mle =
the_est)

  theints <-
boot.ci(the.boot,conf=conf,type=c("norm","basic",
"stud", "perc")) 

  if(logdata) return(c(exp(the_est[1]),exp(theints$norm[2:3])))

  return(c(the_est[1],theints$norm[2:3]))

}

 ```


