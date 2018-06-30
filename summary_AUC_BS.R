rm(list=ls())

library(survival)
library(abind)
source("functions.R")


I = 200 # number of subjects in training dataset
I.test = 100 # number of subjects in testing dataset
n.sim = 2 # number of simulation

# simulated data information
obstime = c(0,3,6,9,12,15, 18,21) # longitudinal measurement time
Tstart = c(6, 9, 12, 15) # landmark time for prediction
deltaT = c(3, 6) # prediction windowns
argvals = obstime/21 # scale the time domain to [0,1]
scenario = "linear_no_missing"

# Table landmark time and predict windows
row.names = cbind(rep(Tstart, each=length(deltaT)), rep(deltaT, length(Tstart)))

models = c("MJM", "MFPCCox")

for(model in models){
  auctrue.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim)
  auc.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim)
  BS.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim)
  MSE.m = matrix(NA, nrow=nrow(row.names), ncol=n.sim) 
  
  file = paste(c("./output/sim_results_", scenario,"_", model), collapse="")
  
  
  for(i.run in 1:n.sim){
    load(paste(c(file, i.run, ".rdata"), collapse=""))
    
    # est K-M curve for BS calculation
    surv = data$surv
    surv.train = surv[surv$id%in%c(1:I), ]
    km = survfit(Surv(time, event)~1, data=surv.train )
    survest = stepfun(km$time, c(1, km$surv))
    
    ith = 0
    for(t in Tstart){

      for(dt in deltaT){
        ith = ith + 1
        tp = t + dt
        
        surv.test = surv[surv$id%in%DP.id[[ith]], ]
        N_vali = nrow(surv.test)
      
        long = data$long
        # The prediction of future multivariate longitudinal outcomes is not available in JMBayes
        tmp.data = long[long$id%in%DP.id[[ith]]&long$obstime<=t, ]
        Y = abind(matrix(tmp.data$Y1, ncol=which(t==obstime), byrow=T), 
                  matrix(tmp.data$Y2, ncol=which(t==obstime), byrow=T),
                  matrix(tmp.data$Y3, ncol=which(t==obstime), byrow=T), along=3)
        erro = abind(matrix(tmp.data$erro[,1], ncol=which(t==obstime), byrow=T), 
                  matrix(tmp.data$erro[,2], ncol=which(t==obstime), byrow=T),
                  matrix(tmp.data$erro[,3], ncol=which(t==obstime), byrow=T), along=3)
        Y = Y-erro  # denoised true longitudinal measures
        
        
        
        #true AUC
        roc = tdROC( X = 1-trueProb[[ith]], Y = timeEvent[[ith]]$time, 
                     delta = timeEvent[[ith]]$event,
                     tau = tp, span = 0.05,
                     nboot = 0, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
        
        auctrue.m[ith, i.run] = roc$AUC$value
        
        #AUC
        roc = tdROC( X = 1-DP.prob[[ith]], Y = timeEvent[[ith]]$time, 
                     delta = timeEvent[[ith]]$event,
                     tau = tp, span = 0.05,
                     nboot = 0, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
        
        auc.m[ith, i.run] = roc$AUC$value
        
        #BS
        D = rep(0, N_vali)
        D[surv.test$time<=tp&surv.test$event==1] = 1
        pi = 1-DP.prob[[ith]]
        
        km_pts = survest(surv.test$time)/survest(t)
        W2 <- D/km_pts
        W1 <- as.numeric(surv.test$time>tp)/(survest(tp)/survest(t))
        W <- W1 + W2
        
        BS_pts <- W * (D - pi)^2
        BS.m[ith, i.run]  = sum(na.omit(BS_pts)) / N_vali
      
        #MSE
        MSE.m[ith, i.run] = mean((DP.long[[ith]][,1:which(t==obstime),] - Y)^2)
        
      }
    }
  }
  
  table = cbind(row.names,  apply(auctrue.m, 1, mean, na.rm=T),  apply(auc.m, 1, mean, na.rm=T), apply(BS.m, 1, mean, na.rm=T) ,
                apply(MSE.m, 1, mean, na.rm=T))
  
  colnames(table) = c("t", "dt", "True.AUC", "AUC", "BS", "MSE")
  
  print(model)
  print(table)
}



