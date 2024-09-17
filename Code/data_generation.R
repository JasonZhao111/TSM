library(PWEALL)
library(survival)


data_generation<-function(n1=200                 #Treatment
                         ,n0=200                 #Control
                         ,rp20=0.75              #Proportion of crossover
){
  #Sample sizes
  n=n1+n0 #Total
  
  tchange=c(0,1,2) #Time-points where the hazards change
  
  rate10=c(0.2,0.2,0.25) #Control hazard before crossover
  rate20=rate10*0.8      #Control hazard after crossover to treatment
  rate30=c(0.4,0.4,0.4)  #Control hazard to crossover
  rate40=rate10*1.5      #Control hazard after crossover and remain in control
  type0=5                #Crossover type=5 (hybrid crossover): rp20 to rate20, (1-rp20) to rate40
  
  
  
  rate11=rate10*0.6      #Treatment hazard before crossover
  rate21=rate11          #Treatment hazard after crossover
  rate31=c(0.2,0.2,0.2)  #Treatment hazard to crossover
  type1=1                #Markov crossover
  
  # Censoring rate
  ratec=c(0.05,0.05,0.05)
  
  #Censoring time data
  tfix=5 #maximum study length
  taur=2 #recruitment window
  u=runif(n)*taur #Randomization time
  tc=rpwe(nr=n,rate=ratec,tchange=tchange)$r
  tc=pmin(tc,tfix-u)
  
  # Compute the (log)hazard ratio
  

  # (log) hazard ratio when crossover is not allowed (true treatment effect)
  log_HR=ovbeta(tfix=tfix,taur=taur,pi1=0.5,
         rate11=rate11,rate21=rate21,rate31=rate31,ratec1=ratec,
         rate10=rate10,rate20=rate20,rate30=rate30,rate40=rate40,ratec0=ratec,
         tchange=tchange,type1=type1,type0=type0,
         rp20=0)$b1
  
  true_HR<-exp(log_HR)
  # Data generation
  
  #Event time data
  abc1=rpwecx(nr=n1,rate1=rate11,rate2=rate21,rate3=rate31,tchange=tchange,type=type1)
  abc0=rpwecx(nr=n0,rate1=rate10,rate2=rate20,rate3=rate30,rate4=rate40,tchange=tchange,type=type0,rp2=rp20)
  
  tpros=c(abc1$rx,abc0$rx) #progression/crossover time
  os=c(abc1$r,abc0$r)      #os time
  pfs=pmin(os,tpros)       #pfs time
  
  
  
  ind.pfs=as.numeric(pfs<=tc) #observed pfs indicator: 1: pfs time is observed, 0: pfs time is censored
  ind.os=as.numeric(os<=tc) #observed os indicator: 1: os time is observed, 0: os time is censored
  cx.to.ind=c(rep(-1,n1),abc0$cxind[,2])   #crossover to treatment indicator: -1: orginally in treatment, 1: crossover to treatment, 0: remain in control
  cx.to.ind[1:n1]=-1
  
  # assemble the dataset
  adata=data.frame(z=c(rep(1,n1),rep(0,n0)), #Treatment indicator
                   pfs=round(pmin(pfs,tc),digits=2), #observed pfs time
                   ind.pfs=ind.pfs,
                   os=round(pmin(os,tc),digits=2), #observed os time
                   ind.os=ind.os,
                   cx.ind=c(abc1$cxind[,1],abc0$cxind[,1])*ind.pfs, #crossover indicator: 1: crossover occured, 0: crossover did not
                   cx.to.ind=cx.to.ind,      #crossover to treatment indicator: -1: orginally in treatment, 1: crossover to treatment, 0: remain in control
                   true_HR= true_HR
                   )
  
  return(adata)
  
}






