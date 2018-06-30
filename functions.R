#####################################################
# functions.R
# Kan Li
# update Jun. 29, 2018
# 
# support functions for data simulation and MFPCCox model
# 
# Functions:
#
# sim_mjm_linear    # simulate Multivariate JM liner longitudinal outcomes
# sim_mjm_nonlinear # simulate Multivariate JM nonliner longitudinal outcomes
# uPACE       # univariate FPCA via PACE
# mfpca.score # MFPC scores
# mfpca.pred  # MFPC longitudinal prediction
# cond.prob   # Risk conditional probability
# functions in tdROC package
# 
#####################################################


# simulation data for multivariate joint model linear 
sim_mjm_linear = function(I, J=8,  obstime = c(0,3,6,9,12,15, 18,21), miss = FALSE, miss.rate = 0.1){
  
  # I : number of subjects
  # J : number of visits
  # obstime: observation times
  # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
  # miss.rate: missing rate.
  
  
  N = I*J
  
  # longitudinal submodel  
  beta0 = c(1.5,2,0.5)
  beta1 = c(2,-1,1)
  betat = c(1.5, -1, 0.6) 
  b.var = c(1,1.5,2)
  e.var = c(1,1,1)
  rho = c(-0.2,0.1,-0.3)
  b.Sigma = diag(b.var)
  b.Sigma[1,2] = b.Sigma[2,1] = sqrt(b.var[1]*b.var[2])*rho[1]
  b.Sigma[1,3] = b.Sigma[3,1] = sqrt(b.var[1]*b.var[3])*rho[2]
  b.Sigma[2,3] = b.Sigma[3,2] = sqrt(b.var[2]*b.var[3])*rho[3]
  
  # sample covariate
  X = rep(rnorm(I, 3, 1), each=J)
  # sample random effect
  ranef = mvrnorm(I, c(0,0,0), b.Sigma)
  id = rep(1:I,each=J)
  ranef = ranef[id,]
  # construct longitudinal submodel
  eta.long = matrix(0, nrow=N, ncol=3)
  for(i in 1:3)
    eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
  
  
  # survival submodel
  gamma1 = -2.5
  alpha = c(0.1, -0.1, 0.2)
  W = rbinom(I, size = 1, prob=.5)
  eta.surv = gamma1*W + c(alpha%*%t(eta.long[!duplicated(id),]))
  
  # simulate survival time
  scale = exp(-7)
  S = runif(I)
  Ti = rep(NA, I)
  alpha.beta = alpha%*%betat
  f = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c(alpha.beta)*t)
    }
    S[i] - exp(-stats::integrate(h, 0, tau)$value)
  }
  f = Vectorize(f)

  for(i in 1:I){
    Ti[i] = uniroot(f, c(0, 100))$root
  }
  
  # simulate true survival probability
  pre.surtime = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c(alpha.beta)*t)
    }
    
    exp(-stats::integrate(h, 0, tau)$value)
  }
  
  true.prob = matrix(NA, nrow=I, ncol=length(obstime[-1]))
  for(i in 1:I){
    ith = 0
    for(tau in obstime[-1]){
      ith = ith + 1
      true.prob[i, ith] = pre.surtime(tau)
    }
  }
  
  colnames(true.prob) = as.character(obstime[-1])
  
  #--------------------------------
  # simulate censor time
  C = runif(I,min=obstime[3], max=obstime[5]+20)
  time <- pmin(Ti, C) #observed time is min of censored and true
  event <- ifelse(time==Ti, 1, 0) #0: censored ; 1: event; 

  # prepare data
  visit = rep(1:J, I)
  obstime = rep(obstime, I) 
  erro = mvrnorm(N, c(0,0,0), diag(e.var))
  Y = matrix(0, nrow=N, ncol=3)
  for(i in 1:3)
    Y[,i] = eta.long[,i] + betat[i]*obstime  + erro[,i]
  
  
  long.all = data.frame(id=id, visit=visit, time = rep(time, each=J), event = rep(event, each=J),
                    Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],obstime=obstime, X=X, ranef=ranef, W  = W[rep(1:I, each=J)], erro=I(erro))
  
  long = long.all
  
  # introduce missing complete at random
  if(miss){
    miss.index = sample(which(long$obstime>obstime[2]), 0.1*N)
    long = long[!c(1:N)%in%miss.index, ]
  }

  surv = data.frame(id = c(1:I),time=time, event=event, W = W, true.prob=I(true.prob))
  
  # remove observations after event or censoring
  long = long[long$obstime<long$time, ]
  
  return(list(long=long, surv=surv, long.all=long.all))
}

# simulation data for multivariate joint model nonlinear 
sim_mjm_nonlinear = function(I, J=8,  obstime = c(0,3,6,9,12,15, 18,21), miss = FALSE, miss.rate = 0.1){
  
  # I : number of subjects
  # J : number of visits
  # obstime: observation times
  # miss: whether introduce missing (missing complete at random) in longitudinal data. Different from drop-out
  # miss.rate: missing rate.
  
  
  N = I*J
  
  
  # longitudinal submodel  
  beta0 = c(1.5,2,0.5)
  beta1 = c(2,-1,1)
  betat = c(1.5, -1, 0.6) 
  b.var = c(1,1.5,2)
  e.var = c(1,1,1)
  rho = c(-0.2,0.1,-0.3)
  b.Sigma = diag(b.var)
  b.Sigma[1,2] = b.Sigma[2,1] = sqrt(b.var[1]*b.var[2])*rho[1]
  b.Sigma[1,3] = b.Sigma[3,1] = sqrt(b.var[1]*b.var[3])*rho[2]
  b.Sigma[2,3] = b.Sigma[3,2] = sqrt(b.var[2]*b.var[3])*rho[3]
  
  # sample covariate
  X = rep(rnorm(I, 3, 1), each=J)
  # sample random effect
  ranef = mvrnorm(I, c(0,0,0), b.Sigma)
  id = rep(1:I,each=J)
  ranef = ranef[id,]
  # construct longitudinal submodel
  eta.long = matrix(0, nrow=N, ncol=3)
  for(i in 1:3)
    eta.long[,i] = beta0[i] + beta1[i]*X + ranef[,i]
  
  
  # survival submodel
  gamma1 = -2.5
  alpha = c(0.1, -0.1, 0.2)
  W = rbinom(I, size = 1, prob=.5)
  eta.surv = gamma1*W + c(alpha%*%t(eta.long[!duplicated(id),]))
  
  # simulate survival time
  scale = exp(-7)
  S = runif(I)
  Ti = rep(NA, I)
  alpha.beta = alpha%*%betat
  c = c(1.2, 0.7, 0.5)
  knots = c(6, 13)
  f = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
                   c[3]*c(alpha.beta)* ifelse(t>=knots[2], t-knots[2], 0) )
    }
    S[i] - exp(-stats::integrate(h, 0, tau)$value)
  }
  f = Vectorize(f)
  
  
  
  for(i in 1:I){
    Ti[i] = uniroot(f, c(0, 100))$root
  }
  
  # simulate true survival probability
  pre.surtime = function(tau){
    h = function(t) {
      scale *exp(eta.surv[i] + c[1]*c(alpha.beta)*ifelse(t<knots[1], t, knots[1]) + c[2]*c(alpha.beta)*ifelse(t<knots[2], pmax(0,(t-knots[1])), knots[2]-knots[1]) +
                   c[3]*c(alpha.beta)* ifelse(t>=knots[2], t-knots[2], 0) ) 
      
    }
    exp(-stats::integrate(h, 0, tau)$value)
  }
  
  true.prob = matrix(NA, nrow=I, ncol=length(obstime[-1]))
  for(i in 1:I){
    ith = 0
    for(tau in obstime[-1]){
      ith = ith + 1
      true.prob[i, ith] = pre.surtime(tau)
    }
  }
  
  colnames(true.prob) = as.character(obstime[-1])
  
  
  #--------------------------------
  # simulate censor time
  C = runif(I,min=obstime[4], max=obstime[5]+25)
  time <- pmin(Ti, C) #observed time is min of censored and true
  event <- ifelse(time==Ti, 1, 0) #0: censored ; 1: event; 
  
  # prepare data
  visit = rep(1:J, I)
  obstime = rep(obstime, I) 
  erro = mvrnorm(N, c(0,0,0), diag(e.var))
  Y = matrix(0, nrow=N, ncol=3)
  for(i in 1:3)
    Y[,i] = eta.long[,i] + betat[i]*(c[1]*ifelse(obstime<knots[1], obstime, knots[1]) + c[2]*ifelse(obstime<knots[2], pmax(0,(obstime-knots[1])), knots[2]-knots[1]) +
                                       c[3]* ifelse(obstime>=knots[2], obstime-knots[2], 0))  + erro[,i]
  
  
  long.all = data.frame(id=id, visit=visit, time = rep(time, each=J), event = rep(event, each=J),
                        Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],obstime=obstime, X=X, ranef=ranef, W  = W[rep(1:I, each=J)], erro=I(erro))
  
  long = long.all
  
  # introduce missing complete at random
  if(miss){
    miss.index = sample(which(long$obstime>obstime[2]), 0.1*N)
    long = long[!c(1:N)%in%miss.index, ]
  }
  
  surv = data.frame(id = c(1:I),time=time, event=event, W = W, true.prob=I(true.prob))
  # remove observations after event or censoring
  long = long[long$obstime<long$time, ]
  
  return(list(long=long, surv=surv, long.all=long.all))
}


# univariate FPCA via principal analysis by conditional estimation(PACE)
uPACE = function(testData, domain, predData=NULL, nbasis = 10, pve = 0.95, npc = NULL){
  
  tmp = funData(domain, testData)
  if(is.null(predData)){
    tmp2 = NULL
  }else{
    tmp2 = funData(domain, predData)
  }
  
  res = PACE(tmp, tmp2, pve=pve, npc= npc, nbasis=nbasis)
  return(res)
} 

# multivariate FPCA based on results from uPACE
mFPCA = function(Xi, phi, p , L){
  
  # eigenanalysis on matrix M
  M = t(Xi)%*%Xi/(I-1)
  eigen.M = eigen(M)
  values = eigen.M$values
  pve = cumsum(values)/sum(values)
  Cms = eigen.M$vectors
  index = unlist(lapply(1:length(L), function(x) rep(x, L[x])))
  
  # MFPCA score
  rho = mfpca.score(Xi, Cms)
  
  # MFPCA eigenfunction
  psis = NULL
  for(j in 1:p){
    psi = NULL
    for(m in 1:dim(Cms)[2]){
      psi = cbind(psi, phi[[j]]%*%Cms[which(index==j),m])
    }
    psis[[j]] = psi
  }
  
  out = list(eigenvalue = values, Cms = Cms, pve = pve, index=index, rho = rho, psis=psis)
  
  return(out)
}

# mfpc score calculation
mfpca.score = function(predXi, Cms){
  rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
  for(i in 1:nrow(predXi)){
    for(m in 1:dim(Cms)[2]){
      rho[i,m] = predXi[i,]%*%Cms[,m]
    }
  }
  return(rho)
}


# mfpc trajectories prediction
mfpca.pred = function(score, meanf, psi, n.rho=NULL){
  p = length(psi)
  n = nrow(score)
  
  if(is.null(n.rho)){
    n.rho = ncol(score)
  }
  
  pred = array(NA, c(n, length(meanf[[1]]), p))
  for(m in 1:p){
    pred[,,m] = matrix(meanf[[m]], nrow=n, ncol =length(meanf[[m]]), byrow = T ) + score[,1:n.rho]%*%t(psi[[m]][, 1:n.rho])
  }
  
  out = pred
  return(out)
}


#risk predict
cond.prob = function(model, newdata, Tstart, Tpred){
  risk.Tstart = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tstart)$surv)
  risk.Tpred = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tpred)$surv)
  return(risk.Tpred/risk.Tstart)
}

################################################################################
##
##                       Functions for R package tdROC
##
################################################################################

calc.AUC <- function( sens, spec ) {
  # Given a vector of sensitivity and specificity, calculate the area under the
  # ROC curve (AUC)
  # Arguments:
  #  -- sens: sensitivity
  #  -- spec: specificity
  # Return:
  #  -- AUC as a numerical scalar
  # NOTE: The AUC is obtained by trapezoidal integration of the area under the
  #       piecewise linear curve obtained by connecting the sensitivity and
  #       specificity.
  o <- order(sens, 1-spec) ;
  y <- sens[o] ;
  x <- 1-spec[o] ;
  x <- c(0, x, 1) ;
  y <- c(0, y, 1) ;
  m <- length(x) ;
  x1 <- x[-m] ;
  x2 <- x[-1] ;
  y1 <- y[-m] ;
  y2 <- y[-1] ;
  AUC <- sum( (y1 + y2)*(x2 - x1)/2 ) ;
  AUC ;
}


calc.kw <- function( X, x0, span ) {
  # Calculate the nearest neighbor kernel weights
  # Arguments:
  #  -- span: the proportion of observations used
  #  -- x0: the x0 around which the kernel weights are calculated
  #  -- X: the vector of all biomarker values in the data
  # Return:
  #  a vector of kernel weights for each element in X
  # NOTE: X must be the vector of ALL biomarker values in the data; it cannot
  #       be any other vector of arbitrary length
  n <- length(X) ;
  tmp0 <- abs( X-x0 ) ;
  tmp1 <- sort( tmp0 ) ;
  tmp2 <- tmp1[ ceiling(n*span) ] ;
  # the cut off that defines the neighborhood
  ans <- as.numeric( tmp0 <= tmp2 ) ;
  ans ;
}


tdROC <- function( X, Y, delta, tau, span, cut.off = NULL,
                   nboot = 0, alpha = 0.05, n.grid = 1000,
                   X.min = NULL, X.max = NULL ) {
  # Calculate the time-dependent sensitivity and specificity and
  # area under the curve (AUC)
  # Arguments:
  #  -- X: the vector of biomarker values
  #  -- Y: time to event
  #  -- delta: indicator of event (1) or censoring (0)
  #  -- tau: the prediction horizon
  #  -- span: the proportion of observations used in calculating kernel weights
  #           (i.e., bandwidth in nearest neighbor method)
  #  -- nboot: number of bootstrap, to be used the variance estimation;
  #            nboot = 0 corresponds to no variance estimation
  #  -- alpha: 1-level of confidence interval, default to be 0.05. It is
  #            used only when nboot > 0
  #  -- n.grid: number of biomarker cut off values used
  #             when calculating the ROC curve
  #  -- cut.off: a vector of biomarker cut.off values at which sensitivity and
  #              specificity will be calculated
  #  -- X.min, X.max: the min and max of X; if NULL, they will be calculated
  #             inside the function; these are not needed for point estimate,
  #             they are needed for bootstrap
  # Return:
  #  -- ROC: a data frame of three columns:
  #             grid, sensitivity and specificity
  #  -- AUC: a data frame of one row and four columns
  #             AUC, standard error of AUC, the lower and upper bootstrap CI
  #  -- prob: a data frame of three columns:
  #             cut.off, sensitivity and specificity
  # NOTE: X, Y and delta must have the same length
  n <- length(X) ;
  positive <- rep(NA, n) ;
  for (i in 1:n) {
    if ( Y[i] > tau ) {
      positive[i] <- 0 ;
    } else {
      if ( delta[i] == 1 ) {
        positive[i] <- 1 ;
      } else {
        kw <- calc.kw( X, X[i], span ) ;
        fm <- survfit( Surv(Y, delta) ~ 1, weights = kw ) ;
        tmp <- summary(fm, times = c(Y[i], tau))$surv ;
        if ( tmp[1] == 0 ) {
          positive[i] <- 1 ;
        } else {
          positive[i] <- 1 - tmp[2]/tmp[1] ;
        }
      }
    }
  }
  negative <- 1 - positive ;
  if ( is.null(X.min) ) { X.min <- min(X) }
  if ( is.null(X.max) ) { X.max <- max(X) }
  grid <- c( -Inf, seq( X.min, X.max, length=n.grid ), Inf ) ;
  sens <- spec <- NULL ;
  for (this.c in grid ) {
    sens <- c( sens, sum(positive*as.numeric(X > this.c))/sum(positive) ) ;
    # sensitivity that incorporates fractional "positive"
    spec <- c( spec, sum(negative*as.numeric(X <= this.c))/sum(negative) ) ;
    # specificity that incorporates fractional "negative"
  }
  ROC <- data.frame( grid = grid,
                     sens = sens,
                     spec = spec ) ;
  AUC <- data.frame( value = calc.AUC( sens, spec ),
                     sd = NA,
                     lower = NA,
                     upper = NA ) ;
  sens <- spec <- NULL ;
  if ( !is.null(cut.off) ) {
    for (this.c in cut.off ) {
      sens <- c( sens, sum(positive*as.numeric(X > this.c))/sum(positive) ) ;
      # sensitivity that incorporates fractional "positive"
      spec <- c( spec, sum(negative*as.numeric(X <= this.c))/sum(negative) ) ;
      # specificity that incorporates fractional "negative"
    }
    prob <- data.frame( cut.off = cut.off,
                        sens = sens,
                        spec = spec ) ;
  } else {
    prob <- NULL ;
  }
  
  if ( nboot > 0 ) {
    # start bootstrap for AUC
    boot.AUC <- rep(NA, nboot) ;
    if ( !is.null(cut.off) ) {
      boot.sens <- matrix( NA, nrow=nboot, ncol=length(cut.off) ) ;
      boot.spec <- matrix( NA, nrow=nboot, ncol=length(cut.off) ) ;
    }
    set.seed(123) ;
    # the random number seed is hardcoded
    for (b in 1:nboot) {
      loc <- sample( x = 1:n, size = n, replace = T ) ;
      X2 <- X[loc] ;
      Y2 <- Y[loc] ;
      delta2 <- delta[loc] ;
      out <- tdROC( X2, Y2, delta2, tau, span, nboot = 0, alpha, n.grid,
                    cut.off = cut.off, X.min = X.min, X.max = X.max ) ;
      boot.AUC[b] <- out$AUC$value ;
      if ( !is.null(cut.off) ) {
        boot.sens[b, ] <- out$prob$sens ;
        boot.spec[b, ] <- out$prob$spec ;
      }
    }
    tmp1 <- sd(boot.AUC) ;
    tmp2 <- as.numeric( quantile( boot.AUC, prob = c(alpha/2, 1-alpha/2) ) ) ;
    AUC$sd <- tmp1 ;
    AUC$lower <- tmp2[1] ;
    AUC$upper <- tmp2[2] ;
    #
    if ( !is.null(cut.off) ) {
      prob$sens.sd <- apply( boot.sens, 2, sd ) ;
      prob$sens.lower <- apply( boot.sens, 2, quantile, prob = alpha/2 ) ;
      prob$sens.upper <- apply( boot.sens, 2, quantile, prob = 1-alpha/2 ) ;
      prob$spec.sd <- apply( boot.spec, 2, sd ) ;
      prob$spec.lower <- apply( boot.spec, 2, quantile, prob = alpha/2 ) ;
      prob$spec.upper <- apply( boot.spec, 2, quantile, prob = 1-alpha/2 ) ;
    } else {
      prob$sens.sd <- NA ;
      prob$sens.lower <- NA ;
      prob$sens.upper <- NA ;
      prob$spec.sd <- NA ;
      prob$spec.lower <- NA ;
      prob$spec.upper <- NA ;
    }
  }
  
  pct.ctrl <- mean( Y > tau ) ;
  pct.case <- mean( Y <= tau & delta == 1 ) ;
  pct.not.sure <- mean( Y <= tau & delta == 0 ) ;
  return( list( ROC = ROC, AUC = AUC, prob = prob,
                pct = c( ctrl = pct.ctrl,
                         case = pct.case,
                         not.sure = pct.not.sure ),
                tau = tau, span = span ) ) ;
}


plot.tdROC <- function( fm ) {
  # Plot the ROC curve estimated by tdROC()
  # Arguments:
  #  -- fm: the object returned by tdROC()
  # Return:
  #  -- a plot of ROC curve
  x <- 1 - fm$ROC$spec ;
  y <- fm$ROC$sens ;
  tmp <- order(x, y) ;
  x2 <- x[tmp] ;
  y2 <- y[tmp] ;
  windows() ;
  plot( x = x2, y = y2, xlab="1-specificity", ylab="sensitivity",
        type="l", lwd = 2, xlim=c(0,1), ylim=c(0,1),
        main=paste("AUC = ", round(fm$AUC$value, 3), " (",
                   round(fm$AUC$lower, 3), ", ",
                   round(fm$AUC$upper, 3), ")", sep="")
  ) ;
  abline(0, 1, col="gray", lwd=2, lty=2) ;
  invisible(0) ;
}


is.monotone <- function(x) {
  # determine if a vector x is monotone (increasing or decreasing) or not
  m <- length(x) ;
  if ( all( x[-m] - x[-1] >= 0 ) | all( x[-m] - x[-1] <= 0 ) ) {
    ans <- TRUE ;
  } else {
    ans <- FALSE ;
  }
  ans ;
}