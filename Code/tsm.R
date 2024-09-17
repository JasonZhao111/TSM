tsm<-function(z,                            #Treatment indicator   1: Treatment, 0: control
              pfs,                          #observed pfs time
              ind.pfs,                      #observed pfs indicator: 1: pfs time is observed, 0: pfs time is censored
              os,                           #observed os time
              ind.os,                       #observed os indicator: 1: os time is observed, 0: os time is censored
              cx.ind,                       #crossover indicator: 1: crossover occured, 0: crossover did not
              cx.to.ind,                    #crossover to treatment indicator: -1: orginally in treatment, 1: crossover to treatment, 0: remain in control
              methods,                      #specify which method, or combination of methods we use
              true_HR="unknown"             #true HR when crossover is not allowed
              )
{
  
  ## load survival package
  library(survival)
  
  ## check for whether or not 100% crossover
  
  pfs_control<-pfs[z==0]                                   # control data pfs
  os_control<-os[z==0]                                     # control data os
  cx.to.ind_control<-cx.to.ind[z==0]                       # 
  
  
  index_progression<- pfs_control<os_control              # control data progression index          
  cp<- sum(index_progression & cx.to.ind_control)/ sum(index_progression)  # crossover proportion
  
  
  ## output matrix
  
  nrow_table<-length(methods)
  
  output_table<-matrix(NA,nrow=nrow_table,ncol=6)
  colnames(output_table)<-c("HR estimation","CI lower","CI upper","standard error","bias","mse")
  
  row.names(output_table)<- rep("name", nrow_table)
  
  k_table=1
  
  
  
  
  ## data preparation
  
  n0<-sum(z==0)                                           # number of observations in control
  n1<-sum(z==1)                                           # number of observations in treatment
  os_treatment<-os[z==1]                                 # treatment data os
  ind.os_control<-ind.os[z==0]                            # control data observed indicator
  ind.os_treatment<-ind.os[z==1]                          # treatment data observed indicator
  ind.pfs_control<-ind.pfs[z==0]                          # control data pfs observed indicator
  group<-z                                                # group index
  
  
  
  # ITT Method
  if ("itt" %in% methods){
    f_itt<-coxph(Surv(time=os,event=ind.os)~z)
    p_itt<-exp(f_itt$coefficients)
    se_itt<-sqrt(f_itt$var)
    lbd_itt<-exp(f_itt$coefficients-1.96*se_itt)
    ubd_itt<-exp(f_itt$coefficients+1.96*se_itt)
    width_itt<-ubd_itt-lbd_itt
    
    output_table[k_table,1]<- p_itt
    output_table[k_table,4]<- se_itt
    output_table[k_table,2]<- lbd_itt
    output_table[k_table,3]<- ubd_itt
    
    if (is.double(true_HR)) {
      bias_itt<- p_itt-true_HR
      square_bias_itt<-bias_itt^2
      mse_itt<-square_bias_itt^2+ f_itt$var
      
      output_table[k_table,5]<- bias_itt
      output_table[k_table,6]<- mse_itt
    }
    
    row.names(output_table)[k_table]<-"itt"
    k_table<-k_table+1
  }
  
  
  ## Censoring at Switch
  
  if ("cas" %in% methods){
    switch_index_cas<- pfs_control < os_control & cx.to.ind_control==1
    
    os_cas<-os_control
    status_cas<-ind.os_control
    
    os_cas[switch_index_cas]<-pmin(os_control, pfs_control)[switch_index_cas]
    status_cas[switch_index_cas]<-pmin(ind.os_control, ind.pfs_control)[switch_index_cas]
    
    
    
    t_cas<-c(os_treatment,os_cas)
    event_cas<-c(ind.os_treatment,status_cas)
    
    
    f_cas<-coxph(Surv(time=t_cas,event=event_cas)~group)
    p_cas<-exp(f_cas$coefficients)
    se_cas<-sqrt(f_cas$var)
    lbd_cas<-exp(f_cas$coefficients-1.96*se_cas)
    ubd_cas<-exp(f_cas$coefficients+1.96*se_cas)
    width_cas<-ubd_cas-lbd_cas
    
    output_table[k_table,1]<- p_cas
    output_table[k_table,4]<- se_cas
    output_table[k_table,2]<- lbd_cas
    output_table[k_table,3]<- ubd_cas
    
    
    if (is.double(true_HR)){
      bias_cas<- p_cas-true_HR
      square_bias_cas<-bias_cas^2
      mse_cas<- f_cas$var+ square_bias_cas
      
      output_table[k_table,5]<- bias_cas
      output_table[k_table,6]<- mse_cas
    }
    
    row.names(output_table)[k_table]<-"cas"
    k_table<-k_table+1
  }
  
  ## exclude at switch 
  if ("eas" %in% methods){
    switch_index_eas<- pfs<os & cx.to.ind==1
    
    t_eas<-os[!switch_index_eas]
    status_eas<-ind.os[!switch_index_eas]
    group_eas<-z[!switch_index_eas]
    
    
    f_eas<-coxph(Surv(time=t_eas,event=status_eas)~group_eas)
    p_eas<-exp(f_eas$coefficients)
    se_eas<-sqrt(f_eas$var)
    
    lbd_eas<-exp(f_eas$coefficients-1.96*se_eas)
    ubd_eas<-exp(f_eas$coefficients+1.96*se_eas)
    width_eas<-ubd_eas-lbd_eas
    
    output_table[k_table,1]<- p_eas
    output_table[k_table,4]<- se_eas
    output_table[k_table,2]<- lbd_eas
    output_table[k_table,3]<- ubd_eas
    
    if (is.double(true_HR)){
      bias_eas<- p_eas-true_HR
      square_bias_eas<-bias_eas^2  
      mse_eas<- f_eas$var+square_bias_eas
      
      output_table[k_table,5]<- bias_eas
      output_table[k_table,6]<- mse_eas
    }
    
    row.names(output_table)[k_table]<-"eas"
    k_table<-k_table+1
  }
  
  # rpsft method
  
  if ("rpsft" %in% methods){
    
    library(rpsftm)
    
    switch_index_control<- pfs_control < os_control & cx.to.ind_control==1
    
    rx_c<-rep(0,n0)
    rx_c[switch_index_control]<-(os_control[switch_index_control]-pfs_control[switch_index_control])/os_control[switch_index_control]
    rx<-c(rep(1,n1),rx_c)
    
    
    f_rpsft<-rpsftm(formula=Surv(os, ind.os) ~ rand(group, rx))
    p_rpsft<-exp(f_rpsft$psi)
    
    
    rpsft_t_control<- (1-rx_c)*os_control + rx_c*p_rpsft*os_control
    rpsft_t<-c(os_treatment,rpsft_t_control)
    
    f_rpsft_hr<-coxph(Surv(time=rpsft_t,event=ind.os)~group)
    
    p_rpsft_hr<-exp(f_rpsft_hr$coefficients)
    
    
    
    if (is.na(p_rpsft)==0 ) {
      
      log_hr_rpsft<-f_rpsft_hr$coefficients
      itt_z <- summary(f_itt)[["coefficients"]][,"z"]
      se_rpsft<-abs(log_hr_rpsft/itt_z)
      p_rpsft_lbd<-exp(log_hr_rpsft-1.96*se_rpsft)
      p_rpsft_ubd<-exp(log_hr_rpsft+1.96*se_rpsft)
      width_rpsft<-p_rpsft_ubd-p_rpsft_lbd
      
      output_table[k_table,1]<- p_rpsft_hr
      output_table[k_table,4]<- se_rpsft
      output_table[k_table,2]<- p_rpsft_lbd
      output_table[k_table,3]<- p_rpsft_ubd
      
      
      if (is.double(true_HR)){
        bias_rpsft <- p_rpsft_hr-true_HR
        square_bias_rpsft<-bias_rpsft^2
        mse_rpsft<-square_bias_rpsft+ se_rpsft^2
        
        output_table[k_table,5]<- bias_rpsft
        output_table[k_table,6]<- mse_rpsft
      }
      
      row.names(output_table)[k_table]<-"rpsft"
      k_table<-k_table+1
    }
  }
  
  
  ## ipcw
  
  if ("ipcw" %in% methods){
    
    switch_index_control<- pfs_control < os_control & cx.to.ind_control==1
    model<-glm(switch_index_control~pfs_control,family="binomial")
    
    ## predicted probability of treatment crossover
    pred_proba <- (1-predict(model,
                             type = "response"))
    
    ## weight cannot be too large
    ipcw_weight<-pmin(c(1/pred_proba[!switch_index_control]),10)
    
    ## ipcw methods remove the crossovered data
    switch_index_eas<- pfs<os & cx.to.ind==1
    t_eas<-os[!switch_index_eas]
    status_eas<-ind.os[!switch_index_eas]
    group_eas<-z[!switch_index_eas]
    
    
    
    f_ipcw<-coxph(Surv(time=t_eas,event=status_eas)~group_eas,weights = c(rep(1,n1),ipcw_weight))
    p_ipcw<-exp(f_ipcw$coefficients)
    se_ipcw<-sqrt(f_ipcw$var)
    
    
    lbd_ipcw<-exp(f_ipcw$coefficients-1.96*se_ipcw)
    ubd_ipcw<-exp(f_ipcw$coefficients+1.96*se_ipcw)
    width_ipcw<-ubd_ipcw-lbd_ipcw
    
    output_table[k_table,1]<- p_ipcw
    output_table[k_table,4]<- se_ipcw
    output_table[k_table,2]<- lbd_ipcw
    output_table[k_table,3]<- ubd_ipcw
    
    
    if(is.double(true_HR)){
      bias_ipcw<- p_ipcw-true_HR
      square_bias_ipcw<-bias_ipcw^2 
      mse_ipcw<-f_ipcw$var+square_bias_ipcw
      
      output_table[k_table,5]<- bias_ipcw
      output_table[k_table,6]<- mse_ipcw
    }
    
    row.names(output_table)[k_table]<-"ipcw"
    k_table<-k_table+1
  }
  
  
  ## treatment as time-dependent variable
  
  if ("ttdv" %in% methods){
    switch_index_control<- pfs_control < os_control & cx.to.ind_control==1
    n_switch<-sum(switch_index_control)
    
    time_start<-c(rep(0,n0),pfs_control[switch_index_control],rep(0,n1))   #starting time
    t_s_c<-os_control                    
    t_s_c[switch_index_control]<-pfs_control[switch_index_control]
    time_stop<-c(t_s_c,os_control[switch_index_control],os_treatment)
    
    s_c<-ind.os_control
    s_c[switch_index_control]<-0
    status_td<-c(s_c, ind.os_control[switch_index_control],ind.os_treatment)
    
    group_td<-c(rep(0,n0),rep(1,n_switch+n1))
    
    
    index_na<-time_start==time_stop
    
    f_td <- coxph(Surv(time_start[!index_na], time_stop[!index_na], status_td[!index_na]) ~ group_td[!index_na])
    p_td<-exp(f_td$coefficients)
    se_td<-sqrt(f_td$var)
    
    lbd_td<-exp(f_td$coefficients-1.96*se_td)
    ubd_td<-exp(f_td$coefficients+1.96*se_td)
    width_td<-ubd_td-lbd_td
    
    output_table[k_table,1]<- p_td
    output_table[k_table,4]<- se_td
    output_table[k_table,2]<- lbd_td
    output_table[k_table,3]<- ubd_td
    
    
    if (is.double(true_HR)) {
      bias_td<- p_td-true_HR
      square_bias_td<-bias_td^2
      mse_td<-f_td$var+square_bias_td
      
      output_table[k_table,5]<- bias_td
      output_table[k_table,6]<- mse_td
    }
    
    row.names(output_table)[k_table]<-"ttdv"
    k_table<-k_table+1
  }
  
  
  
  
  
  
  ## two-stage method
  
  if ("tse" %in% methods){
    
    if (cp !=1) {
      switch_index_control<- pfs_control < os_control & cx.to.ind_control==1
      os_progression<- os_control[index_progression]
      pfs_progression<- pfs_control[index_progression]
      
      ind.os_progression<- ind.os_control[index_progression]
      cx.to.ind_progression<- cx.to.ind_control[index_progression]
      
      t_p<-os_progression-pfs_progression
      
      
      aft_fit<-survreg(Surv(t_p ,ind.os_progression) ~ cx.to.ind_progression, dist = "weibull")
      
      exp_phi<-exp(aft_fit$coefficients[2])
      
      # adjust for survival time
      
      os_adjusted<-os_control
      
      os_adjusted[switch_index_control]<-(os_control[switch_index_control]-pfs_control[switch_index_control])/exp_phi+pfs_control[switch_index_control]
      
      
      
      t_tse<-c(os_treatment,os_adjusted)
      
      
      f_tse<-coxph(Surv(time=t_tse,event=ind.os)~z)
      p_tse<-exp(f_tse$coefficients)
      
      
      
      ## bootstrap to calculate variance
      
      bootstrap_iter=200
      p_bootstrap<-rep(0,200)
      
      
      
      # bootstrap sample
      
      for (i in 1:bootstrap_iter) {
        adata<-data.frame(cbind(z,pfs,ind.pfs,ind.os,os,cx.to.ind))
        index_bootstrap<-sample(1:nrow(adata),nrow(adata),replace = TRUE)
        adata_bootstrap<-adata[index_bootstrap,]
        
        bootstrap_control_data<-adata_bootstrap[adata_bootstrap$z==0,]
        adata_treatment<-adata_bootstrap[adata_bootstrap$z==1,]
        
        n_control<-nrow(bootstrap_control_data)
        n_treatment<-nrow(adata_treatment)
        
        switch_index_control<- (bootstrap_control_data$pfs<bootstrap_control_data$os) & (bootstrap_control_data$cx.to.ind==1)
        
        
        
        #index_bootstrap<-sample(1:n_control, n_control, replace = TRUE)
        #bootstrap_control_data<-adata_control[index_bootstrap,]
        progression_adata_control<- bootstrap_control_data[bootstrap_control_data$pfs<bootstrap_control_data$os,]
        
        t_p<-progression_adata_control$os-progression_adata_control$pfs
        
        
        aft_fit<-survreg(Surv(t_p ,ind.os) ~ cx.to.ind,
                         data = progression_adata_control, dist = "weibull")
        
        exp_phi<-exp(aft_fit$coefficients[2])
        
        # adjust for survival time
        
        adata_control_adjusted<-bootstrap_control_data
        
        adata_control_adjusted$os[switch_index_control]<-(bootstrap_control_data$os[switch_index_control]-bootstrap_control_data$pfs[switch_index_control])/exp_phi+bootstrap_control_data$pfs[switch_index_control]
        
        #index_treatment<-sample(1:n_control, n_treatment, replace = TRUE)
        t_e<-adata_treatment$os
        status_e<-adata_treatment$ind.os
        
        
        t_tse<-c(t_e,adata_control_adjusted$os)
        
        status_tse<-c(status_e,adata_control_adjusted$ind.os)
        
        group_bootstrap<- c(rep(1,n_treatment),rep(0,n_control))
        f_tse<-coxph(Surv(time=t_tse,event=status_tse)~group_bootstrap)
        p_tse_bootstrap<-exp(f_tse$coefficients)
        p_bootstrap[i]<-p_tse_bootstrap
      }
      
      
      
    }
    
    
    
    if (cp==1 ) {
      
     
  
      switch_index_control<- pfs_control < os_control & cx.to.ind_control==1
      os_progression<- os_control[index_progression]
      pfs_progression<- pfs_control[index_progression]
      
      ind.os_progression<- ind.os_control[index_progression]
      cx.to.ind_progression<- cx.to.ind_control[index_progression]
      
      t_p<-os_progression-pfs_progression
      
      adata<-data.frame(cbind(z,pfs,ind.pfs,ind.os,os,cx.to.ind))
      adata_control<- adata[z==0,]
      adata_control_adjusted<-adata[z==0,]
  

  
      phi=0 # first time iteration, no adjustment at the first time
      phi_aft=1 # make the when loop work at the first iteration
  
      tol_aft<- 1e-4
  
     iter_tse=1
      while (abs(phi-phi_aft)>tol_aft) {
    
    if (iter_tse>1) phi<-phi_aft
    exp_phi<-exp(phi)
    
    adata_control_adjusted$os[switch_index_control]<-(adata_control$os[switch_index_control]-adata_control$pfs[switch_index_control])/exp_phi+adata_control$pfs[switch_index_control]
    
    t_tse<-c(os_treatment,adata_control_adjusted$os)

    ## remove the os=0, weibull dist has to be >0
    index_model<-t_tse!=0
    aft_tse<-survreg(Surv(time=t_tse[index_model],event=ind.os[index_model])~group[index_model])
    phi_aft<-aft_tse$coefficients[2]
    iter_tse<-iter_tse+1
  }

 
  
  adata_control_adjusted$os[switch_index_control]<-(adata_control$os[switch_index_control]-adata_control$pfs[switch_index_control])/exp_phi+adata_control$pfs[switch_index_control]
  t_tse<-c(os_treatment,adata_control_adjusted$os)

  f_tse<-coxph(Surv(time=t_tse,event=ind.os)~group)
  p_tse<-exp(f_tse$coefficients)
  
  
  ## bootstrap to calculate variance
  
  bootstrap_iter=200
  p_bootstrap<-rep(0,200)
  
  
  
  # bootstrap sample
  
  for (i in 1:bootstrap_iter) {
    index_bootstrap<-sample(1:nrow(adata),nrow(adata),replace = TRUE)
    adata_bootstrap<-adata[index_bootstrap,]
    
    bootstrap_control_data<-adata_bootstrap[adata_bootstrap$z==0,]
    adata_treatment<-adata_bootstrap[adata_bootstrap$z==1,]
    
    n_control<-nrow(bootstrap_control_data)
    n_treatment<-nrow(adata_treatment)
    
    group_bootstrap<-c(rep(1,n_treatment),rep(0,n_control))
    switch_index_control<- (bootstrap_control_data$pfs<bootstrap_control_data$os) & (bootstrap_control_data$cx.to.ind==1)
    
    
    
    progression_adata_control<- bootstrap_control_data[bootstrap_control_data$pfs<bootstrap_control_data$os,]
    
    t_p<-progression_adata_control$os-progression_adata_control$pfs
    
    adata_control_adjusted<-bootstrap_control_data
    
    phi=0 # first time iteration, no adjustment at the first time
    phi_aft=1 # make the when loop work at the first iteration
    
    tol_aft<- 1e-4
    
    iter_tse=1
    while (abs(phi-phi_aft)>tol_aft) {
      
      if (iter_tse>1) phi<-phi_aft 
      exp_phi<-exp(phi)
      
      adata_control_adjusted$os[switch_index_control]<-(bootstrap_control_data$os[switch_index_control]-bootstrap_control_data$pfs[switch_index_control])/exp_phi+bootstrap_control_data$pfs[switch_index_control]
      
      t_e<-adata_treatment$os
      status_e<-adata_treatment$ind.os
      
      
      t_tse<-c(t_e,adata_control_adjusted$os)
      status_tse<-c(status_e,adata_control_adjusted$ind.os)
      
      
      
      
      ## remove the os=0, weibull dist has to be >0
      index_model<-t_tse!=0
      
      aft_tse<-survreg(Surv(time=t_tse[index_model],event=status_tse[index_model])~group_bootstrap[index_model])
      phi_aft<-aft_tse$coefficients[2]
      iter_tse<-iter_tse+1
    }
    
    # adjust for survival time
    
    adata_control_adjusted<-bootstrap_control_data
    
    adata_control_adjusted$os[switch_index_control]<-(bootstrap_control_data$os[switch_index_control]-bootstrap_control_data$pfs[switch_index_control])/exp_phi+bootstrap_control_data$pfs[switch_index_control]
    
    #index_treatment<-sample(1:n_control, n_treatment, replace = TRUE)
    t_e<-adata_treatment$os
    status_e<-adata_treatment$ind.os
    
    
    t_tse<-c(t_e,adata_control_adjusted$os)
    status_tse<-c(status_e,adata_control_adjusted$ind.os)
    
    
    f_tse<-coxph(Surv(time=t_tse,event=status_tse)~group_bootstrap)
    p_tse_bootstrap<-exp(f_tse$coefficients)
    p_bootstrap[i]<-p_tse_bootstrap
  }
  
    }
    
     se_tse<- sd(p_bootstrap)
     lbd_tse<-quantile(p_bootstrap,0.025)
     ubd_tse<-quantile(p_bootstrap,0.975)
     width_tse<-ubd_tse-lbd_tse


     output_table[k_table,1]<- p_tse
     output_table[k_table,4]<- se_tse
     output_table[k_table,2]<- lbd_tse
     output_table[k_table,3]<- ubd_tse


    if (is.double(true_HR)){
        bias_tse <- p_tse-true_HR
        square_bias_tse<-bias_tse^2
        mse_tse<-square_bias_tse+ var(p_bootstrap)
  
        output_table[k_table,5]<- bias_tse
        output_table[k_table,6]<- mse_tse
    }
    
    
    row.names(output_table)[k_table]<-"tse"
    k_table<-k_table+1
  
  }  
  
  
  ## our three-state model / BIMM method.: two scenarios cp==1 or cp==0
  
  if ("bimm" %in% methods) {
    
    library(rstan)

    
    cutpoints<-seq(from=0,to=max(os_control),length=6)
    switch_index_control<- pfs_control < os_control & cx.to.ind_control==1
    
    # specify how many CPU cores needed (correspond to how many chains)
    
    n_chain=2             #speify number of chain
    options(mc.cores=n_chain)    
    
    
    # check for 100% crossover
    
    if (cp!=1) {
      dat<-list(K=length(cutpoints)-1,
                cutpoints=cutpoints,
                N=n0,
                cens_t=1-ind.os_control,
                cens_tx=1-as.numeric(pfs_control<os_control),
                crossover_type=2-cx.to.ind_control,
                t=os_control,
                tx=pfs_control
      )
      

      
      fit=stan(file = 'ph_estimable.stan', data = dat,chains = n_chain)
      
      
      burns=1001
      l=2000
      hr_3s<-rep(0,n_chain*(l-burns+1))
      model_var<-rep(0,n_chain*(l-burns+1))
      
      
      
      ## extract lambdas
      
      for (chain in 1:n_chain) {
        fixeffect <- fit@sim$samples[[chain]][grep("lambda", names(fit@sim$samples[[chain]]))]
        t_regenerate_matrix<-matrix(0,(l-burns+1),n0)
    
        for (tt in burns:l) {
          lam<-matrix(0,4,length(cutpoints)-1)
          
          ## store lambdas estimate
          
          #lam[1,] is estimated lambda1
          #lam[2,] is estimated lambda2*
          #lam[3,] is estimated lambda2
          #lam[4,] is estimated lambda3  
          k=1  
          for (i in 1:4) {
            for (j in 1:(length(cutpoints)-1)) {
              lam[i,j]<-fixeffect[[k]][tt]
              k<-k+1
            }
          }
          
          #lam[2,] is estimated lambda2*
          #lam[3,] is estimated lambda2
          lam2_ns<-lam[3,]
          lam2_s<-lam[2,]
          
          t_regenerate<-os_control
          
          for (ii in 1:n0) {
            
            if (switch_index_control[ii]==1){
              
              t_s<-os_control[ii]-pfs_control[ii]
              s_prob<-st(t_s,lambda = lam2_s,cutpoints = cutpoints[1:(length(cutpoints)-1)])
              
              t_regenerate[ii]<-pfs_control[ii]+uniroot(st,interval = c(0,2e09),lambda=lam2_ns,cutpoints=cutpoints[1:(length(cutpoints)-1)],root_finding=s_prob)$root
              
            }
            
          }  
          t_regenerate_matrix[tt-burns+1,]<-t_regenerate
        }
        
        
        for (jj in 1: (l-burns+1)) {
          t_regenerate_all<-c(os_treatment,t_regenerate_matrix[jj,])
          status_regenerate<-ind.os
          f_3s<-coxph(Surv(time=t_regenerate_all,event=ind.os)~group)
          p_3s<-exp(f_3s$coefficients)
          hr_3s[(l-burns+1)*(chain-1)+jj]=p_3s
          model_var[(l-burns+1)*(chain-1)+jj]=f_3s$var
        }
      } 
    }
     
    if (cp==1) {
      
      ## assuming hazard ratio before and after progression are the same
      pfs.fit<- coxph(Surv(time=pfs,event=ind.pfs)~z)
      p_pfs<-exp(pfs.fit$coefficients)
      
      
      dat<-list(K=length(cutpoints)-1,
                cutpoints=cutpoints,
                N=n0,
                cens_t=1- ind.os_control,
                cens_tx=1-as.numeric(pfs_control<os_control),
                t=os_control,
                tx=pfs_control
      )
      
      fit=stan(file = 'pem3.stan', data = dat,chains = n_chain)
      
      l=2000
      burns=1501
      
      hr_3s<-rep(0,n_chain*(l-burns+1))
      model_var<-rep(0,n_chain*(l-burns+1))
      
      length_each_chain= l-burns+1
      ## extract lambdas
      
      for (chain in 1:2) {
        fixeffect <- fit@sim$samples[[chain]][grep("lambda", names(fit@sim$samples[[chain]]))]
        t_regenerate_matrix<-matrix(0,(l-burns+1),n0)
        
        
        
        for (tt in burns:l) {
          lam<-matrix(0,3,length(cutpoints)-1)
          
          ## store lambdas estimate
          
          #lam[1,] is estimated lambda1
          #lam[2,] is estimated lambda2*
          #lam[3,] is estimated lambda3  
          k=1  
          for (i in 1:3) {
            for (j in 1:(length(cutpoints)-1)) {
              lam[i,j]<-fixeffect[[k]][tt]
              k<-k+1
            }
          }
          
          itertime=1
          beta_star<-p_pfs
          
          while (itertime<=100) {
            
            itertime=itertime+1
            

            
            #lam[2,] is estimated lambda2*
            #lam[3,] is estimated lambda2
            lam2_s<-lam[2,]
            lam2_ns<-lam[2,]/beta_star
            
            t_regenerate<-os_control
            
            for (ii in 1:n0) {
              
              if (switch_index_control[ii]==1){
                
                t_s<-os_control[ii]-pfs_control[ii]
                s_prob<-st(t_s,lambda = lam2_s,cutpoints = cutpoints[1:(length(cutpoints)-1)])
                

                t_regenerate[ii]<-pfs_control[ii]+uniroot(st,interval = c(0,2e09),lambda=lam2_ns,cutpoints=cutpoints[1:(length(cutpoints)-1)],root_finding=s_prob)$root
              }
              
            }  
            
            t_regenerate_all<-c(os_treatment,t_regenerate)
            f_3s<-coxph(Surv(time=t_regenerate_all,event=ind.os)~z)
            p_3s<-exp(f_3s$coefficients)
            beta<- p_3s
            
            ## stop when 
            if (abs(beta-beta_star)<1e-6) itertime=105
            beta_star=beta
          }
          
          hr_3s[length_each_chain*(chain-1)+tt-burns+1]<-beta_star
          model_var[length_each_chain*(chain-1)+tt-burns+1]<-f_3s$var
        }
        
      }
    }
    
    
    p_3s<- mean(hr_3s)
    se_3s<-sqrt(mean(model_var)+var(log(hr_3s)))
    lbd_threestate=exp(mean(log(hr_3s))-1.96*se_3s)
    ubd_threestate=exp(mean(log(hr_3s))+1.96*se_3s)
    
    
    
    output_table[k_table,1]<- p_3s
    output_table[k_table,4]<- se_3s
    output_table[k_table,2]<- lbd_threestate
    output_table[k_table,3]<- ubd_threestate
    
    
    if (is.double(true_HR)){
      bias_3s <- p_3s-true_HR
      square_bias_3s<-bias_3s^2
      mse_3s<-square_bias_3s+ se_3s^2
      
      output_table[k_table,5]<- bias_3s
      output_table[k_table,6]<- mse_3s
    }
    
    row.names(output_table)[k_table]<-"bimm"
  }      
  
  
  
  # plotting
  
  library(miscTools)
  
  credplot.gg <- function(d){
    # d is a data frame with 4 columns
    # d$x gives variable names
    # d$y gives center point
    # d$ylo gives lower limits
    # d$yhi gives upper limits
    require(ggplot2)
    p <- ggplot(d, aes(x = factor(x, level = level_order), y=y, ymin=ylo, ymax=yhi))+
      geom_pointrange()+
      coord_flip()+
      xlab('')+
      ylab("Hazard Ratio")
    return(p)
  }
  
  
  row.names(output_table)<-toupper(row.names(output_table))
  method<-row.names(output_table)
  hr_est<-output_table[,1]
  low<-output_table[,2]
  up<-output_table[,3]
  
  level_order <- method[length(method):1]       ## setting the order
  
  
  d<-data.frame(x=method,y=hr_est,ylo=low,yhi=up)
  
  p<-credplot.gg(d)
  
  if (is.double(true_HR)){
            print(p +
            theme(axis.text=element_text(size=24),
             axis.title=element_text(size=24,face="bold")))
    
    
  }
  
  if (is.double(true_HR)){
    print(p + geom_hline(yintercept = true_HR,col="blue",lty=2)+theme(axis.text=element_text(size=24),
                                                                axis.title=element_text(size=24,face="bold")))
    
    
  }

  return(output_table)
}



## Helper functions for the BIMM method

st<-function(t,lambda,cutpoints,root_finding=0){
  
  
  ## number of pieces
  k<-length(lambda);
  
  ## H is cumulative hazard
  
  H<-0
  
  
  ## length in each pieces, last piece to infinity
  
  length<-rep(0,k-1)
  
  
  
  for (j in 1:(k-1)) {
    length[j]<-cutpoints[j+1]-cutpoints[j]
  }
  
  
  for (i in 1:(k-1)) {
    if (t >= cutpoints[i+1]) {
      H<-H+lambda[i]*length[i]
    } else if (cutpoints[i] <= t && t <cutpoints[i+1]) {
      H<-H+lambda[i]*(t-cutpoints[i])
    } 
  }
  
  if (t>cutpoints[k]){
    H<-H+lambda[k]*(t-cutpoints[k])
  }
  
  ## return the survival prob, root_finding is used to generate variable
  
  s<-exp(-H)-root_finding
  return(s)
}


## generate function

pw_simulation<-function(n,lambda,cutpoints){
  
  time<-rep(0,n)
  
  F<-runif(n)
  
  for (i in 1:n) {
    time[i]<-uniroot(st,interval = c(0,1e09),lambda=lambda,cutpoints=cutpoints,root_finding=F[i])$root
  }
  return(time)
}




  












