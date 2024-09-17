source("data_generation.R")
source("tsm.R")

## Generate data
set.seed(123456)
adata_sample <- data_generation(n1=200,n0=200,rp20=0.75)
attach(adata_sample)
head(adata_sample)

# true treatment effect
true_HR <- true_HR[1]


## Example 1: true HR known, choose different subset of methods

# TSM methods to estimate the treatment effect
methods=c("itt","cas","eas","rpsft","ipcw","tse","ttdv")

tsm_example_1 <- tsm(z=z,                            #Treatment indicator   1: Treatment, 0: control
                     pfs=pfs,                        #observed pfs time
                     ind.pfs=ind.pfs,                #observed pfs indicator: 1: pfs time is observed, 0: pfs time is censored
                     os=os,                          #observed os time
                     ind.os=ind.os,                  #observed os indicator: 1: os time is observed, 0: os time is censored
                     cx.ind=cx.ind,                  #crossover indicator: 1: crossover occurred, 0: crossover did not
                     cx.to.ind=cx.to.ind,            #crossover to treatment indicator: -1: originally in treatment, 1: crossover to treatment, 0: remain in control
                     methods=methods,                #specify which method, or combination of methods we use
                     true_HR=true_HR                 #true HR when crossover is not allowed
)

# Results
tsm_example_1

# To extract different methods' point estimates for treatment effect
tsm_example_1[,1]
# Estimates and 95% uncertainty intervals
tsm_example_1[,c(1:3)]

# To extract results for a specific method, e.g. RPSFT
tsm_example_1[c("RPSFT"), ]
# e.g. results for TSE
tsm_example_1[c("TSE"), ]


## Example 2: true HR unknown, 100 percent crossover

## Generate data
set.seed(123456)
adata_sample <- data_generation(n1=200,n0=200,rp20=1)
attach(adata_sample)

# TSM methods to estimate the treaatment effect
methods=c("itt","cas","eas","rpsft","ipcw","tse","ttdv")

tsm_example_2 <- tsm(z=z,                            #Treatment indicator   1: Treatment, 0: control
                     pfs=pfs,                        #observed pfs time
                     ind.pfs=ind.pfs,                #observed pfs indicator: 1: pfs time is observed, 0: pfs time is censored
                     os=os,                          #observed os time
                     ind.os=ind.os,                  #observed os indicator: 1: os time is observed, 0: os time is censored
                     cx.ind=cx.ind,                  #crossover indicator: 1: crossover occurred, 0: crossover did not
                     cx.to.ind=cx.to.ind,            #crossover to treatment indicator: -1: originally in treatment, 1: crossover to treatment, 0: remain in control
                     methods=methods,                #specify which method, or combination of methods we use
                     true_HR="unknown"               #true HR when crossover is not allowed. Either a number of "unknown"
)

# Since true_HR="unknown", the bias and mse columns show "NA"
tsm_example_2

# Methods' point estimates for treatment effect
tsm_example_2[,1]
# Estimates and 95% uncertainty intervals
tsm_example_2[,c(1:3)]

# To extract results for a specific method, e.g. TTDV
tsm_example_2[c("TTDV"), ]


## Example 3: try BIMM and TSE with 75% crossover

# Generate data
set.seed(123456)
adata_sample <- data_generation(n1=200,n0=200,rp20=0.75)
attach(adata_sample)
true_HR<-true_HR[1]

# ITT, TSE, and BIMM methods
methods<-c("itt","tse","bimm")

tsm_example_3 <- tsm(z=z,                            #Treatment indicator   1: Treatment, 0: control
                     pfs=pfs,                        #observed pfs time
                     ind.pfs=ind.pfs,                #observed pfs indicator: 1: pfs time is observed, 0: pfs time is censored
                     os=os,                          #observed os time
                     ind.os=ind.os,                  #observed os indicator: 1: os time is observed, 0: os time is censored
                     cx.ind=cx.ind,                  #crossover indicator: 1: crossover occurred, 0: crossover did not
                     cx.to.ind=cx.to.ind,            #crossover to treatment indicator: -1: originally in treatment, 1: crossover to treatment, 0: remain in control
                     methods=methods,                #specify which method, or combination of methods we use
                     true_HR=true_HR                 #true HR when crossover is not allowed
)

# Results
tsm_example_3

# Methods' point estimates for treatment effect
tsm_example_3[,1]
# Estimates and 95% uncertainty intervals
tsm_example_3[,c(1:3)]

# To extract results for a specific method, e.g. BIMM
tsm_example_3[c("BIMM"), ]
