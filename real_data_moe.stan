// Stan model for tensile strength given observed knots

data{
  
  int<lower = 1> N;       // total number of lumbers
  int<lower = 2> J;       // number of partition
  int<lower = 1> Kmax;     // maximum number of knots 
  
  vector<lower = 0>[N] y_obs;        // the observed data vector 
  int min_index[N];                     // index of the observed minimum 
  
  real<lower = 0> D_array[J, Kmax, N];  // stacked D matrices
  real<lower = 0> Z_array[Kmax, 1, N]; // stacked Z matrices
  real<lower = -1, upper = 1> Edge_array[Kmax, 1, N]; // stacked Edge Indicator 
  
  real<lower= 0> dmax; // influential distance
  int<lower = 0> K_vec[N]; // number of knots for each specimen
  real<lower=0> MOE_vec[N]; // vector of the observed moe 

} 


parameters{
  // parameters for the AR process and weighted effects.
  // AR(1) parameters
  real eta0;
  real eta1;
  real<lower = 0, upper = 1> rho;
  real<lower = 0> sigma;
  // knot-effects parameters
  real<lower = 0> beta;
  real<lower = 0> gamma0; // scale parameter for regular knots
  real<lower = 0> gamma1; // scale parameter for edge knots
  
  // latent y values
  matrix<lower = 0> [J-1, N] y_lat_raw; 
  
}

transformed parameters{
  
  //  real mu = alpha/(1-rho);
  
  real margin_sigma =  sigma/sqrt(1-rho^2);
  
  matrix<lower=0> [J-1, N] y_lat;
  
  for(i in 1:N){
    for(j in 1:(J-1)){
      y_lat[j,i] = y_lat_raw[j,i] + y_obs[i];
    }
  }
  
}

model{
  // prior
  eta0 ~ normal(0,10);
  eta1 ~ normal(0,10);
  beta ~ normal(0, 1);
  rho ~ normal(.5,.5);
  sigma ~ cauchy(0,5);
  gamma0 ~ normal(0,1);
  gamma1 ~ normal(0,1);
  
  // likelihood 
  for (i in 1:N){
    int K_i; 
    real MOE_i;
    vector[J] weighted_effects;       // weighted_effects for the current lumber
    K_i = K_vec[i];                   // number of knots in that piece of lumber 
    MOE_i = MOE_vec[i];               // MOE for lumber i 
    
                                      // get the weighted effects                        
    for(j in 1:J){                    // loop over the cell j 
      real weighted_effects_j;        // weighted_effects_j - for the current cell j
      weighted_effects_j = 0;     
      for(k in 1:K_i){                // loop over the K_i knots 
        real distance_k;
        real edge_k;
        distance_k = D_array[j,k,i];  // lumber i, cell j, knots k
        edge_k = Edge_array[k,1,i];   // Edge indicator for lumber i, knots  
        if(distance_k <= dmax){       // check if the knot k is within influential distance
          if(edge_k==0){              // edge_k ==0 -> not edge knot
            weighted_effects_j += gamma0*exp(-beta*distance_k) * Z_array[k,1,i];
          } else if (edge_k == 1){    // edge_k ==1 -> is edge knot
            weighted_effects_j += gamma1*exp(-beta*distance_k) * Z_array[k,1,i];
          }
        }
      }
      weighted_effects[j] = weighted_effects_j;  // weighted_effects: WE_1, .... WE_J
    }                      
    if(min_index[i]==1){
      y_obs[i] ~ normal((eta0+eta1*MOE_i) - weighted_effects[1], margin_sigma);
    } else {
      y_lat[1,i] ~ normal((eta0+eta1*MOE_i) - weighted_effects[1], margin_sigma);
    }
    
    for(j in 2:J){
      if(min_index[i] == j){
        y_obs[i] ~  normal((1-rho)*(eta0+eta1*MOE_i) + rho * y_lat[j-1,i] - weighted_effects[j] + rho*weighted_effects[j-1], sigma); 
      } 
      else{
        if(min_index[i]== j-1){
          y_lat[j-1,i] ~ normal((1-rho)*(eta0+eta1*MOE_i) + rho * y_obs[i] - weighted_effects[j] + rho*weighted_effects[j-1], sigma); 
        } else if(min_index[i] < j-1){
          y_lat[j-1,i] ~ normal((1-rho)*(eta0+eta1*MOE_i) + rho * y_lat[j-2,i] - weighted_effects[j] + rho*weighted_effects[j-1], sigma); 
        } else {
          y_lat[j,i] ~ normal((1-rho)*(eta0+eta1*MOE_i) + rho * y_lat[j-1,i] - weighted_effects[j] + rho*weighted_effects[j-1], sigma); 
        }
      }
    }
  }
  
}














