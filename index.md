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

get\_est = function(est\_rt,est\_os,sd\_rt,sd\_os,B=NA,conf=.95){

  var\_os = sd\_os^2

  var\_rt = sd\_rt^2

  bias = est\_os-est\_rt

  w =  (bias^2+var\_os)/(bias^2+var\_os+var\_rt)

  est = w\*est\_rt+(1-w)\*est\_os

  if(is.na(B)) return(est)

  # parametric bootstrap

  data = c(est\_rt=est\_rt,est\_os=est\_os,sd\_rt=sd\_rt,sd\_os=sd\_os)

  est = function(data){

    est\_rt = data[1]

    est\_os = data[2]

    var\_rt = data[3]^2

    var\_os = data[4]^2

    bias = est\_os-est\_rt

    w =  (bias^2+var\_os)/(bias^2+var\_os+var\_rt)

    return(w\*est\_rt+(1-w)\*est\_os)

  }

  the\_est = est(data)

  random\_data = function(data,mle=the\_est){

    est\_par = mle[1]

    est\_bias = data[2]-data[1]

    sd\_rt = data[3]

    sd\_os = data[4]

    rt\_b = rnorm(1,mean=est\_par,sd=sd\_rt)

    os\_b = rnorm(1,mean=est\_par+est\_bias,sd=sd\_os)

    return(c(rt\_b,os\_b,sd\_rt,sd\_os))

  }

  the.boot = boot(data, est, R = B, sim = &quot;parametric&quot;,

                   ran.gen = random\_data, mle = mean(the\_est))

  theints = boot.ci(the.boot,conf=conf,type=c(&quot;norm&quot;,&quot;basic&quot;, &quot;stud&quot;, &quot;perc&quot;))

  # assuming 0 is true parameters do intervals cover 0?

  return(c(the\_est[1],theints$norm[2:3]))

}

##  est\_cons - function to calculate basic estimate given in Section 2.2  of manuscript tapering the weight for the trial toward 1.

##  parameters:

##  est\_rt: treatment effect estimate from randomized trial

##  est\_os: treatment effect estimate from observational study

##  sd\_rt: standard error of treatment effect estimate from randomized trial

##  st\_os: standard error of treatment effect estimate from observational study

##  B = number bootstrap resamples

est\_cons = function(est\_rt,est\_os,sd\_rt,sd\_os,B=NA,conf=.95){

  if(is.na(B)){

    var\_os = sd\_os^2

    var\_rt = sd\_rt^2

    bias = est\_os-est\_rt

    bias = max(abs(bias-2\*sqrt(var\_rt+var\_os)),abs(bias+2\*sqrt(var\_rt+var\_os)))

    w =  (bias^2+var\_os)/(bias^2+var\_os+var\_rt)

    return(w\*est\_rt+(1-w)\*est\_os)

  }

  data = c(est\_rt=est\_rt,est\_os=est\_os,sd\_rt=sd\_rt,sd\_os=sd\_os)

  est = function(data){

    est\_rt = data[1]

    est\_os = data[2]

    var\_rt = data[3]^2

    var\_os = data[4]^2

    bias = est\_os-est\_rt

    bias = max(abs(bias-2\*sqrt(var\_rt+var\_os)),abs(bias+2\*sqrt(var\_rt+var\_os)))

    w =  (bias^2+var\_os)/(bias^2+var\_os+var\_rt)

    return(w\*est\_rt+(1-w)\*est\_os)

  }

  the\_est = est(data)

  random\_data = function(data,mle=the\_est){

    est\_par = mle[1]

    est\_bias = data[2]-data[1]

    sd\_rt = data[3]

    sd\_os = data[4]

    rt\_b = rnorm(1,mean=est\_par,sd=sd\_rt)

    os\_b = rnorm(1,mean=est\_par+est\_bias,sd=sd\_os)

    return(c(rt\_b,os\_b,sd\_rt,sd\_os))

  }

  the.boot = boot(data, est, R = B, sim = &quot;parametric&quot;,

                   ran.gen = random\_data, mle = mean(the\_est))

  theints = boot.ci(the.boot,conf=conf,type=c(&quot;norm&quot;,&quot;basic&quot;, &quot;stud&quot;, &quot;perc&quot;))

  # assuming 0 is true parameters do intervals cover 0?

  return(c(the\_est[1],theints$norm[2:3]))

}

#  fixed effects meta analysis (theta\_vec,se\_vec are vectors)

fe\_meta = function(theta\_vec,se\_vec){

  weights = 1/se\_vec^2/(sum(1/se\_vec^2))

  return(c(sum(weights\*theta\_vec),sum(weights^2\*se\_vec^2)))

}

#  random effects meta analysis using DerSimmion/Laird estimate for between study variance (theta\_vec,se\_vec are vectors)



re\_meta = function(theta\_vec,se\_vec){

  F = fe\_meta(theta\_vec,se\_vec)[1]

  Q = sum(se\_vec^-2\*(theta\_vec-F)^2)

  sigmab2 = max(0, (Q-(length(theta\_vec)-1))/(sum(se\_vec^-2)-sum(se\_vec^-4)/sum(se\_vec^-2)))

  weights = 1/(se\_vec^2+sigmab2)/(sum(1/(se\_vec^2+sigmab2)))

  return(c(sum(weights\*theta\_vec),sum(weights^2\*(se\_vec^2+sigmab2))))

}



#  Fixed-bias meta analysis (from 3.1).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  This function allows negative weights for the observational studies

##  est\_rt: treatment effect vector from collection of randomized trials

##  est\_os: treatment effect vector from collection of observational studies

##  sd\_rt: standard error of treatment effects from randomized trials

##  st\_os: standard error of treatment effects from observational studies

##  B = number bootstrap resamples

##  logdata=TRUE (if the est\_rt and est\_os are given on log-scale (as might be true for logistic regression).  If logdata=TRUE, effects are exponentiated before outputing

est\_fixedeffect = function(est\_rt,est\_os,sd\_rt,sd\_os,B=1000,logdata=TRUE,conf=.95){

  if(is.na(B)){

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat

    var\_os = sd\_os^2

    var\_rt = sd\_rt^2

    Dmat = 1\*(c(rep(0,length(est\_rt)),est\_bias)%\*%t(c(rep(0,length(est\_rt)),est\_bias))+diag(c(sd\_rt^2,sd\_os^2)))

    bvec = c(1)

    meq = 1

    dvec = c(rep(0,length(est\_rt)),rep(0,length(est\_os)))

    Amat = t(matrix(c(rep(1,length(est\_rt)),rep(1,length(est\_os))),nrow=length(bvec)))

    weights = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w\_rt = weights[1:length(est\_rt)]

    w\_os = weights[(length(est\_rt)+1):length(weights)]

    if(logdata) return(exp((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt))))

    return((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt)))

  }

  N\_rt = length(est\_rt)

  N\_os = length(est\_os)

  data = c(est\_rt,sd\_rt,est\_os,sd\_os,N\_rt,N\_os)

  get\_est = function(data){

    N\_os = data[length(data)]

    N\_rt = data[length(data)-1]

    est\_rt = data[1:N\_rt]

    est\_os = data[(2\*N\_rt+1):(2\*N\_rt+N\_os)]

    sd\_rt = data[(N\_rt+1):(2\*N\_rt)]

    sd\_os = data[(2\*N\_rt+N\_os+1):(2\*N\_rt+2\*N\_os)]

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat

    var\_os = sd\_os^2

    var\_rt = sd\_rt^2

    Dmat = 1\*(c(rep(0,length(est\_rt)),est\_bias)%\*%t(c(rep(0,length(est\_rt)),est\_bias))+diag(c(sd\_rt^2,sd\_os^2)))

    bvec = c(1)

    meq = 1

    dvec = c(rep(0,length(est\_rt)),rep(0,length(est\_os)))

    Amat = t(matrix(c(rep(1,length(est\_rt)),rep(1,length(est\_os))),nrow=length(bvec)))

    weights = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w\_rt = weights[1:length(est\_rt)]

    w\_os = weights[(length(est\_rt)+1):length(weights)]

    return((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt)))

  }

  the\_est = get\_est(data)

  random\_data = function(data,mle=the\_est){



    N\_os = data[length(data)]

    N\_rt = data[length(data)-1]

    est\_rt = data[1:N\_rt]

    est\_os = data[(2\*N\_rt+1):(2\*N\_rt+N\_os)]

    sd\_rt = data[(N\_rt+1):(2\*N\_rt)]

    sd\_os = data[(2\*N\_rt+N\_os+1):(2\*N\_rt+2\*N\_os)]

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat

    rt\_b = rnorm(length(sd\_rt),mean=mle,sd=sd\_rt)

    os\_b = rnorm(length(sd\_os),mean=mle+est\_bias,sd=sd\_os)

    return(c(rt\_b,sd\_rt,os\_b,sd\_os,N\_rt,N\_os))

  }

  the.boot = boot(data, get\_est, R = B, sim = &quot;parametric&quot;,

                   ran.gen = random\_data, mle = the\_est)

  theints = boot.ci(the.boot,conf=conf,type=c(&quot;norm&quot;,&quot;basic&quot;, &quot;stud&quot;, &quot;perc&quot;))

  if(logdata) return(c(exp(the\_est),exp(theints$norm[2:3])))

  return(c(the\_est,theints$norm[2:3]))

}

#  Fixed-bias meta analysis (from 3.1).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  This function does not allow negative weights for the observational studies

##  est\_rt: treatment effect vector from collection of randomized trials

##  est\_os: treatment effect vector from collection of observational studies

##  sd\_rt: standard error of treatment effects from randomized trials

##  st\_os: standard error of treatment effects from observational studies

##  B = number bootstrap resamples

##  logdata=TRUE (if the est\_rt and est\_os are given on log-scale (as might be true for logistic regression).  If logdata=TRUE, effects are exponentiated before outputing



est\_fixedeffect\_nn = function(est\_rt,est\_os,sd\_rt,sd\_os,B=1000,logdata=TRUE,conf=.95){

  if(is.na(B)){

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat

    var\_os = sd\_os^2

    var\_rt = sd\_rt^2

    Dmat = 1\*(c(rep(0,length(est\_rt)),est\_bias)%\*%t(c(rep(0,length(est\_rt)),est\_bias))+diag(c(sd\_rt^2,sd\_os^2)))

    bvec = c(1,rep(0,length(est\_rt)),rep(0,length(est\_os)))

    meq = 1

    dvec = c(rep(0,length(est\_rt)),rep(0,length(est\_os)))

    Amat = t(matrix(rbind( c(rep(1,length(est\_rt)),rep(1,length(est\_os))),diag(1,length(est\_rt)+length(est\_os))),nrow=length(bvec)))

    weights = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w\_rt = weights[1:length(est\_rt)]

    w\_os = weights[(length(est\_rt)+1):length(weights)]

    if(logdata) return(exp((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt))))

    return((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt)))

  }

  N\_rt = length(est\_rt)

  N\_os = length(est\_os)

  data = c(est\_rt,sd\_rt,est\_os,sd\_os,N\_rt,N\_os)

  get\_est = function(data){

    N\_os = data[length(data)]

    N\_rt = data[length(data)-1]

    est\_rt = data[1:N\_rt]

    est\_os = data[(2\*N\_rt+1):(2\*N\_rt+N\_os)]

    sd\_rt = data[(N\_rt+1):(2\*N\_rt)]

    sd\_os = data[(2\*N\_rt+N\_os+1):(2\*N\_rt+2\*N\_os)]

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat

    var\_os = sd\_os^2

    var\_rt = sd\_rt^2

    Dmat = 1\*(c(rep(0,length(est\_rt)),est\_bias)%\*%t(c(rep(0,length(est\_rt)),est\_bias))+diag(c(sd\_rt^2,sd\_os^2)))

    bvec = c(1,rep(0,length(est\_rt)),rep(0,length(est\_os)))

    meq = 1

    dvec = c(rep(0,length(est\_rt)),rep(0,length(est\_os)))

    Amat = t(matrix(rbind( c(rep(1,length(est\_rt)),rep(1,length(est\_os))),diag(1,length(est\_rt)+length(est\_os))),nrow=length(bvec)))

    weights = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w\_rt = weights[1:length(est\_rt)]

    w\_os = weights[(length(est\_rt)+1):length(weights)]

    return((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt)))

  }

  the\_est = get\_est(data)

  random\_data = function(data,mle=the\_est){



    N\_os = data[length(data)]

    N\_rt = data[length(data)-1]

    est\_rt = data[1:N\_rt]

    est\_os = data[(2\*N\_rt+1):(2\*N\_rt+N\_os)]

    sd\_rt = data[(N\_rt+1):(2\*N\_rt)]

    sd\_os = data[(2\*N\_rt+N\_os+1):(2\*N\_rt+2\*N\_os)]

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat

    rt\_b = rnorm(length(sd\_rt),mean=mle,sd=sd\_rt)

    os\_b = rnorm(length(sd\_os),mean=mle+est\_bias,sd=sd\_os)

    return(c(rt\_b,sd\_rt,os\_b,sd\_os,N\_rt,N\_os))

  }

  the.boot = boot(data, get\_est, R = B, sim = &quot;parametric&quot;,

                   ran.gen = random\_data, mle = the\_est)

  theints = boot.ci(the.boot,conf=conf,type=c(&quot;norm&quot;,&quot;basic&quot;, &quot;stud&quot;, &quot;perc&quot;))

  if(logdata) return(c(exp(the\_est),exp(theints$norm[2:3])))

  return(c(the\_est,theints$norm[2:3]))

}

#  Random-bias meta analysis (from 3.1) using DerSimmion/Laird estimate for between study variance (theta\_vec,se\_vec are vectors).  The estimate for the randomized trials is calculated using a fixed effects meta analysis.  This function does not allow negative weights for the trials

##  est\_rt: treatment effect vector from collection of randomized trials

##  est\_os: treatment effect vector from collection of observational studies

##  sd\_rt: standard error of treatment effects from randomized trials

##  st\_os: standard error of treatment effects from observational studies

##  B = number bootstrap resamples

##  logdata=TRUE (if the est\_rt and est\_os are given on log-scale (as might be true for logistic regression).  If logdata=TRUE, effects are exponentiated before outputing

est\_randomeffect = function(est\_rt,est\_os,sd\_rt,sd\_os,B=1000,logdata=TRUE,conf=.95){

  if(is.na(B)){

    N\_os = length(est\_os)

    N\_rt = length(est\_rt)

    out = fe\_meta(est\_rt,sd\_rt)

    theta\_hat = out[1]

    theta\_hat\_os =  fe\_meta(est\_os,sd\_os)[1]

    Q = sum(sd\_os^-2\*(est\_os-theta\_hat\_os)^2)

    sigmab2 = max(0, (Q-(length(est\_os)-1))/(sum(sd\_os^-2)-sum(sd\_os^-4)/sum(sd\_os^-2)))

    est\_bias = rep(fe\_meta(est\_os,sqrt(sd\_os^2+sigmab2))[1] - theta\_hat,length(sd\_os))

    var\_os = sd\_os^2 + sigmab2

    var\_rt = sd\_rt^2

    Dmat = 1\*(c(rep(0,length(est\_rt)),est\_bias)%\*%t(c(rep(0,length(est\_rt)),est\_bias))+diag(c(sd\_rt^2,var\_os)))

    bvec = c(1)

    meq = 1

    dvec = c(rep(0,length(est\_rt)),rep(0,length(est\_os)))

    Amat = t(matrix(c(rep(1,length(est\_rt)),rep(1,length(est\_os))),nrow=length(bvec)))

    weights = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w\_rt = weights[1:length(est\_rt)]

    w\_os = weights[(length(est\_rt)+1):length(weights)]

    if(logdata) return(exp((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt))))

    return((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt)))

  }

  N\_rt = length(est\_rt)

  N\_os = length(est\_os)

  data = c(est\_rt,sd\_rt,est\_os,sd\_os,N\_rt,N\_os)

  est = function(data){

    N\_os = data[length(data)]

    N\_rt = data[length(data)-1]

    est\_rt = data[1:N\_rt]

    est\_os = data[(2\*N\_rt+1):(2\*N\_rt+N\_os)]

    sd\_rt = data[(N\_rt+1):(2\*N\_rt)]

    sd\_os = data[(2\*N\_rt+N\_os+1):(2\*N\_rt+2\*N\_os)]

    out = fe\_meta(est\_rt,sd\_rt)

    theta\_hat = out[1]

    theta\_hat\_os =  fe\_meta(est\_os,sd\_os)[1]

    Q = sum(sd\_os^-2\*(est\_os-theta\_hat\_os)^2)

    sigmab2 = max(0, (Q-(length(est\_os)-1))/(sum(sd\_os^-2)-sum(sd\_os^-4)/sum(sd\_os^-2)))

    est\_bias = rep(fe\_meta(est\_os,sqrt(sd\_os^2+sigmab2))[1] - theta\_hat,length(sd\_os))

    var\_os = sd\_os^2 + sigmab2

    var\_rt = sd\_rt^2

    Dmat = 1\*(c(rep(0,length(est\_rt)),est\_bias)%\*%t(c(rep(0,length(est\_rt)),est\_bias))+diag(c(sd\_rt^2,var\_os)))

    bvec = c(1)

    meq = 1

    dvec = c(rep(0,length(est\_rt)),rep(0,length(est\_os)))

    Amat = t(matrix(c(rep(1,length(est\_rt)),rep(1,length(est\_os))),nrow=length(bvec)))

    weights = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE)$solution

    w\_rt = weights[1:length(est\_rt)]

    w\_os = weights[(length(est\_rt)+1):length(weights)]

    return(c((sum(w\_rt\*est\_rt)+sum(w\_os\*est\_os))/(sum(w\_os)+sum(w\_rt)),sigmab2 ))

  }

  the\_est = est(data)

  random\_data = function(data,mle=the\_est){

    N\_os = data[length(data)]

    N\_rt = data[length(data)-1]

    est\_rt = data[1:N\_rt]

    est\_os = data[(2\*N\_rt+1):(2\*N\_rt+N\_os)]

    sd\_rt = data[(N\_rt+1):(2\*N\_rt)]

    sd\_os = data[(2\*N\_rt+N\_os+1):(2\*N\_rt+2\*N\_os)]

    theta\_hat = fe\_meta(est\_rt,sd\_rt)[1]

    est\_bias = est\_os-theta\_hat+rnorm(N\_os,mean=0,sd=sqrt(mle[2]))

    rt\_b = rnorm(length(sd\_rt),mean=mle[1],sd=sd\_rt)

    os\_b = rnorm(length(sd\_os),mean=mle[1]+est\_bias,sd=sd\_os)

    return(c(rt\_b,sd\_rt,os\_b,sd\_os,N\_rt,N\_os))

  }

  the.boot = boot(data, est, R = B, sim = &quot;parametric&quot;,

                   ran.gen = random\_data, mle = the\_est)

  theints = boot.ci(the.boot,conf=conf,type=c(&quot;norm&quot;,&quot;basic&quot;, &quot;stud&quot;, &quot;perc&quot;))

  if(logdata) return(c(exp(the\_est[1]),exp(theints$norm[2:3])))

  return(c(the\_est[1],theints$norm[2:3]))

}

##  function to convert a 95% confidence interval (on OR scale) to the corresponding standard error on the log-OR scale

getse = function(l1,l2){

  log1 = log(l1)

  log2 = log(l2)

  a = ((log1+log2)/2)

  sel = (log2-log1)/(2\*1.96)

  return(c(a,sel))

}
