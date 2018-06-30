#####################################################
# sim_MFPCCox.R
# Kan Li
# update Jun. 29, 2018
# 
# simulation study to illustrate the implementation 
# of MJM mdoel
#
#####################################################


rm(list=ls())

source("functions.R")
library(JMbayes)
library(MASS)
library(abind)

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

for(i.run in 1:n.sim){
  print(i.run)
  set.seed(123+i.run)
  
  # simulate data
  data = sim_mjm_linear((I+I.test), miss = F)  # three longitudinal outcomes, linear trajectories, without missing
  #data = sim_mjm_nonlinear((I+I.test), miss = FALSE) # three longitudinal outcomes, nonlinear trajectories, without missing
  long = data$long  # longitudinal data
  surv = data$surv  # survival data
  
  # split data to training and testing
  long.train = long[long$id%in%c(1:I), ]
  long.test = long[!long$id%in%c(1:I), ]
  surv.train = surv[surv$id%in%c(1:I), ]
  surv.test = surv[!surv$id%in%c(1:I), ]
  
  # Fit longitudinal submodels
  MixedModelFit = mvglmer(list(Y1 ~ X + obstime + (1 | id),
                                Y2 ~ X + obstime + (1 | id), 
                                Y3 ~ X + obstime + (1 | id)), data =long.train,
                           families = list(gaussian, gaussian, gaussian))
  
  # Fit survival submodel
  CoxFit = coxph(Surv(time, event) ~ W , data = surv.train, model = TRUE)
  
  # Fit MJM
  JMFit = mvJointModelBayes(MixedModelFit, CoxFit, timeVar = "obstime")
  postMeans = JMFit$statistics
  
  DP.id = DP.prob = DP.long = timeEvent = trueProb = NULL
  ith = 0
  # prediction for different time windowes
  for(t in Tstart){
    tmp.id = surv.test[surv.test$time>t, "id"]
    tmp.surv.data = surv.test[surv.test$time>t, ]
    tmp.data = long.test[long.test$id%in%tmp.id&long.test$obstime<=t, ]
    for(dt in deltaT){
      ith = ith + 1
      DP.id[[ith]] = tmp.id
      timeEvent[[ith]] = tmp.surv.data[, c("time", "event")] 
      tmp.pred = survfitJM(JMFit, tmp.data, survTimes= (t+dt))
      trueProb [[ith]] = tmp.surv.data$true.prob[, (which((t+dt)==obstime)-1)]
      DP.prob[[ith]] = do.call(rbind, lapply(1:length(tmp.id), function(x) tmp.pred[[1]][[x]]))[,2]
      DP.long[[ith]] = abind(do.call(rbind, lapply(1:length(tmp.id), function(x) tmp.pred$fitted.y[[x]][[1]])),
                        do.call(rbind, lapply(1:length(tmp.id), function(x) tmp.pred$fitted.y[[x]][[2]])),
                        do.call(rbind, lapply(1:length(tmp.id), function(x) tmp.pred$fitted.y[[x]][[3]])), along=3)
    }
  }
  
  save(data,  trueProb, postMeans, DP.id,DP.prob, DP.long, timeEvent, file=paste(c("./output/sim_results_", scenario, "_MJM", i.run, ".rdata"), collapse=""))

}

print(proc.time() - ptm)

# n.sim = 2
# > print(proc.time() - ptm)
# user  system elapsed 
# 65.00    1.03  401.77
