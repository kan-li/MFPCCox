#####################################################
# sim_MFPCCox.R
# Kan Li
# update Jun. 29, 2018
# 
# simulation study to illustrate the implementation 
# of MFPCCox mdoel
#
#####################################################

rm(list=ls())

source("functions.R")
library(MASS)
library(refund)
library(MFPCA)


I = 200 # number of subjects in training dataset
I.test = 100 # number of subjects in testing dataset
n.sim = 2 # number of simulation

# simulated data information
obstime = c(0,3,6,9,12,15, 18,21) # longitudinal measurement time
Tstart = c(6, 9, 12, 15) # landmark time for prediction
deltaT = c(3, 6) # prediction windowns
argvals = obstime/21 # scale the time domain to [0,1]
scenario = "linear_no_missing"


ptm = proc.time()

# main part for simulation
for(i.run in 1:n.sim){
  print(i.run)
  set.seed(123+i.run)
  
  # simulate data
  data = sim_mjm_linear((I+I.test), miss = F)  # three longitudinal outcomes, linear trajectories, without missing
  #data = sim_mjm_nonlinear((I+I.test), miss = F) # three longitudinal outcomes, nonlinear trajectories, without missing
  long = data$long  # longitudinal data
  surv = data$surv  # survival data
  
  # split data to training and testing
  long.train = long[long$id%in%c(1:I), ]
  long.test = long[!long$id%in%c(1:I), ]
  surv.train = surv[surv$id%in%c(1:I), ]
  surv.test = surv[!surv$id%in%c(1:I), ]
  
  # subject ids
  patID = surv$id
  nPat = length(patID)
  
  # transfer longitudinal outcomes from long to wide
  multivar = array(NA, c((I+I.test), length(obstime), 3))
  for(i in 1:nPat){
    visits = which(obstime %in% (long$obstime[long$id == patID[i]]))
    multivar[i,visits, 1] = long$Y1[long$id == patID[i]]
    multivar[i,visits, 2] = long$Y2[long$id == patID[i]]
    multivar[i,visits, 3]= long$Y3[long$id == patID[i]]
  }
  
  multivar.train = multivar[1:I, , ]
  
  # univariate FPCA via PACE
  Xi.train = L = phi.train = meanFun.train =  NULL
  
  for(p in 1:3){
    tmp.ufpca = uPACE(multivar.train[,,p], argvals, nbasis=3)
    Xi.train = cbind(Xi.train, tmp.ufpca$scores) # FPC scores
    L = c(L, dim(tmp.ufpca$scores)[2])
    phi.train[[p]] = t(tmp.ufpca$functions@X) # FPC eigenfunctions
    meanFun.train[[p]] = tmp.ufpca$mu@X # estimated mean functions
  }
  
  # multivariate FPCA
  mFPCA.train = mFPCA(Xi=Xi.train, phi=phi.train, p=3, L=L )
  rho.train = mFPCA.train$rho  #MFPC scores
  pve = mFPCA.train$pve
  psi = mFPCA.train$psi
  Cms = mFPCA.train$Cms
  
  # survival model 
  surv.train$rho = I(rho.train)
  CoxFit = coxph(Surv(time, event) ~ W + rho , data = surv.train, model = TRUE)
  
  
  # dynamic prediction
  DP.id = DP.prob = DP.long = timeEvent= trueProb = NULL
  ith = 0
  for(t in Tstart){
    tmp.id = surv.test[surv.test$time>t, "id"]  # subjects still event-free at landmark time 
    tmp.surv.data = surv.test[surv.test$time>t, ] # filter the data
    tmp.data = multivar[tmp.id, , ] # subset longitudinal outcomes
    tmp.data[,-c(1:which(t==obstime)),] = NA  # retain the measurements before landmark time
    
    # univariate FPC 
    Xi.test = NULL
    fit.test = array(NA, dim(tmp.data))
    for(p in 1:3){
      tmp.ufpca = uPACE(multivar.train[,,p], argvals, tmp.data[,,p], nbasis=3)
      Xi.test = cbind(Xi.test, tmp.ufpca$scores) # dynamic FPC scores for test subjects 
    }
    
    # estimate MFPC scores for test subjects
    rho.test = mfpca.score(Xi.test, Cms)
    tmp.surv.data$rho = rho.test
    
    # predict longitudinal trajectories 
    long.pred = mfpca.pred(rho.test, meanFun.train, psi)
    
    
    # risk prediction for different time windowes
    for(dt in deltaT){
      ith = ith + 1
      DP.id[[ith]] = tmp.id
      DP.long[[ith]] = long.pred
      timeEvent[[ith]] = tmp.surv.data[, c("time", "event")] # true event time and even indicator
      trueProb [[ith]] = tmp.surv.data$true.prob[, (which((t+dt)==obstime)-1)] # true risk 
      DP.prob[[ith]] = as.numeric(summary(survfit(CoxFit, newdata = tmp.surv.data, se.fit = F, conf.int = F), times = (t+dt))$surv)
    }
  }
  
  save(data, trueProb, Xi.train, rho.train, pve, rho.test, DP.id, DP.prob, DP.long, timeEvent, file=paste(c("./output/sim_results_", scenario,"_MFPCCox", i.run, ".rdata"), collapse=""))

}


print(proc.time() - ptm)

# n.sim = 2
# > print(proc.time() - ptm)
# user  system elapsed 
# 1.78    0.03    2.29