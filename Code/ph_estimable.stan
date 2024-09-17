
 data {

     // Number of pieces
     int<lower=0> K;
     // Cutopoints on time
     //  cutpoints[1] = 0
     //  max(event time) < cutpoints[K+1] < Inf
     //  intervals of cutpoints are closed
     real cutpoints[K+1];
     
     // number of observations
     int<lower=0> N;
     
     // censor index, 1 denotes censor.
     int<lower=0,upper=1> cens_t[N];
     
     // crossover index, 1 denotes crossvoer. 0 denote not crossover
     int<lower=0,upper=1> cens_tx[N];
     
     // crossover type, take value 1 and 2. 1 correspond to lambda_21, 2 correspond to lambda_22. (\lambda_2 and \lambda_2*)
     int crossover_type[N];
     
     // observed event time 
     real t[N];
     
     // observed crossover time ( tx>t is meaningless )
     real tx[N];
 }
 
 transformed data {
    
    // length in each group may differ 
    real length[K];
    
    // time after crossover
    real t2[N];
    
    for (j in 1:K){
       length[j]=cutpoints[j+1]-cutpoints[j];
    }
    
    for (i in 1:N){
       t2[i]=t[i]-tx[i];
    }
    
    
 }
 
 parameters {
     // Baseline hazards, \lambda_2 is estimable
     real<lower=0> lambda1[K];
     real<lower=0> lambda21[K];
     real<lower=0> lambda22[K];
     real<lower=0> lambda3[K];
 }
 
 transformed parameters {

 }
 
 model {
     // prior setting
     
     for (i in 1:K){
         lambda1[K]~ gamma(1,2);
         lambda21[K]~gamma(1,2);
         lambda22[K]~gamma(1,2);
         lambda3[K]~gamma(1,2);
         
     }

     
     
     // compute the log of complete likelihood
     // Loop over pieces of time
     for (i in 1:N) {
       
         // k = 1,2,...,K
         // cutpoints[1] = 0
         // cutpoints[K+1] > max event time
         // Likelihood contribution
         
         if(cens_tx[i]==0){

             // contribute for logs1(tx)+logs3(tx)+log lambda3(tx) part
            for (k in 1:K) {
             // Everyone will contribute to the survival part.
             if (tx[i] >= cutpoints[k+1]) {
                 // If surviving beyond the end of the interval,
                 // contribute survival throughout the interval.
                 target += -lambda1[k] * length[k];
                 target += -lambda3[k] * length[k];
                 //
             } else if (cutpoints[k] <= tx[i] && tx[i] < cutpoints[k+1]) {
                 // If ending follow up during the interval,
                 // contribute survival until the end of follow up.
                 target += -lambda1[k] * (tx[i] - cutpoints[k]);
                 target += -lambda3[k] * (tx[i] - cutpoints[k]);
                 //
                 // Event individuals also contribute to the hazard part.
                target += log(lambda3[k]);
               } 
         }
         
           // contribute for logs2(t-tx)
             if (crossover_type[i]==1){
                for (k in 1:K) {
             // Everyone will contribute to the survival part.
             if ( t2[i]>= cutpoints[k+1]) {
                 // If surviving beyond the end of the interval,
                 // contribute survival throughout the interval.
                 target += -lambda21[k] * length[k];
                 //
             } else if (cutpoints[k] <= t2[i] && t2[i] < cutpoints[k+1]) {
                 // If ending follow up during the interval,
                 // contribute survival until the end of follow up.
                 target += -lambda21[k] * (t2[i] - cutpoints[k]);
                 //
                 // Event individuals also contribute to the hazard part.
                 if (cens_t[i]==0){
                        target += log(lambda21[k]); 
                 }
               } 
         }
                   
             }
             
            if (crossover_type[i]==2){
                for (k in 1:K) {
             // Everyone will contribute to the survival part.
             if ( t2[i]>= cutpoints[k+1]) {
                 // If surviving beyond the end of the interval,
                 // contribute survival throughout the interval.
                 target += -lambda22[k] * length[k];
                 //
             } else if (cutpoints[k] <= t2[i] && t2[i] < cutpoints[k+1]) {
                 // If ending follow up during the interval,
                 // contribute survival until the end of follow up.
                 target += -lambda22[k] * (t2[i] - cutpoints[k]);
                 //
                 // Event individuals also contribute to the hazard part.
                 if (cens_t[i]==0){
                        target += log(lambda22[k]); 
                 }
               } 
         }
                   
             }
         
         
         
         }
         
         
         // log s1(t)+ log s3(t) for everyone not switched
         if (cens_tx[i] ==1){
                    for (k in 1:K) {
             // Everyone will contribute to the survival part.
             if (t[i] >= cutpoints[k+1]) {
                 // If surviving beyond the end of the interval,
                 // contribute survival throughout the interval.
                 target += -lambda1[k] * length[k];
                 target += -lambda3[k] * length[k];
                 //
             } else if (cutpoints[k] <= t[i] && t[i] < cutpoints[k+1]) {
                 // If ending follow up during the interval,
                 // contribute survival until the end of follow up.
                 target += -lambda1[k] * (t[i] - cutpoints[k]);
                 target += -lambda3[k] * (t[i] - cutpoints[k]);
                 //
                 // Event individuals also contribute to the hazard part.
                 if (cens_t[i] == 0) {
                     target += log(lambda1[k]);
                 }
             } 
         }
           
         }

     }
 }


