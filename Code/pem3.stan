
 data {

     // Number of pieces
     int<lower=0> K;
     // Cutopoints on time
     //  cutpoints[1] = 0
     //  max(event time) < cutpoints[K+1] < Inf
     //  K+1 elements
     real cutpoints[K+1];
     //
     int<lower=0> N;
     
     // censor index for t
     int<lower=0,upper=1> cens_t[N];
     
     // censor index for tx
     int<lower=0,upper=1> cens_tx[N];
     
     real t[N];
     real tx[N];
 }
 
 transformed data {
    
    // length in each group may differ 
    real length[K];
    
    real t2[N];
    
    for (j in 1:K){
       length[j]=cutpoints[j+1]-cutpoints[j];
    }
    
    for (i in 1:N){
       t2[i]=t[i]-tx[i];
    }
    
    
 }
 
 parameters {
     // Baseline hazards
     real<lower=0> lambda1[K];
     real<lower=0> lambda2[K];
     real<lower=0> lambda3[K];
 }
 
 transformed parameters {

    
 }
 
 model {
     // prior setting
     
     for (i in 1:K){
         lambda1[K]~ gamma(1,2);
         lambda2[K]~ gamma(1,2);
         lambda3[K]~gamma(1,2);
         
     }

     
     
     
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
           
             for (k in 1:K) {
             // Everyone will contribute to the survival part.
             if ( t2[i]>= cutpoints[k+1]) {
                 // If surviving beyond the end of the interval,
                 // contribute survival throughout the interval.
                 target += -lambda2[k] * length[k];
                 //
             } else if (cutpoints[k] <= t2[i] && t2[i] < cutpoints[k+1]) {
                 // If ending follow up during the interval,
                 // contribute survival until the end of follow up.
                 target += -lambda2[k] * (t2[i] - cutpoints[k]);
                 //
                 // Event individuals also contribute to the hazard part.
                 if (cens_t[i]==0){
                        target += log(lambda2[k]); 
                 }
               } 
         }
         
         
         
         }
         
         
         // log s1(t)+ log s3(t) for everyone with tx censored
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


