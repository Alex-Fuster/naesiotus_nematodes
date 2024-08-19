data {
  int N_spp;
  array[N_spp] int<lower=1, upper=N_spp> sp_id;
  
  // habitat data
  int N_habitat;
  array[N_habitat] int habitat_arboreal;   // Total counts of arid habitat for each species
  array[N_habitat] int total_hab;  //total counts of habitat (arboreal or terrestrial)
  array[N_habitat] int sp_index_hab; 
  
  // vegetation zone data
  int N_veg;
  array[N_veg] int veg_arid;
  array[N_veg] int total_veg;
  array[N_veg] int sp_index_veg;
  
}
parameters {
  real mu_arboreal;
  real<lower=0> sd_arboreal;
  vector[N_spp] arboreal_prob_logit;
  
  real mu_arid;
  real<lower=0> sd_arid;
  vector[N_spp] arid_prob_logit;
}
transformed parameters {
  vector[N_spp] arboreal_prob = inv_logit(arboreal_prob_logit);
  vector[N_spp] arid_prob = inv_logit(arid_prob_logit);
}
model {
  // Hierarchical priors for the prob of microhabitat
  mu_arboreal ~ std_normal();
  sd_arboreal ~ exponential(1);
  arboreal_prob_logit ~ normal(mu_arboreal, sd_arboreal);
  // Likelihood for habitat observations
  habitat_arboreal ~ binomial_logit(total_hab, arboreal_prob_logit[sp_index_hab]);
  
  // Hierarchical priors for vegetation zone
  mu_arid ~ std_normal();
  sd_arid ~ exponential(1);
  arid_prob_logit ~ normal(mu_arid, sd_arid);
  // Likelihood for vegetation zone
  veg_arid ~ binomial_logit(total_veg, arid_prob_logit[sp_index_veg]);
  
  
}
