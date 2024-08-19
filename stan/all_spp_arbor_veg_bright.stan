data {
  // Species and Island data
  int N_spp;
  int N_island;
  array[N_spp] int<lower=1, upper=N_spp> sp_id;
  vector[N_spp] island_area;
  vector[N_spp] island_age;
  array[N_spp] int<lower=1,upper=N_island> island_index_spp;
  
  // habitat data
  int N_habitat;
  array[N_habitat] int habitat_arboreal;   // Total counts of arid habitat for each species
  array[N_habitat] int total_hab;  //total counts of habitat (arboreal or terrestrial)
  array[N_habitat] int<lower=1,upper=N_spp> sp_index_hab; 
  
  // vegetation zone data
  int N_veg;
  array[N_veg] int veg_arid;
  array[N_veg] int total_veg;
  array[N_veg] int<lower=1,upper=N_spp> sp_index_veg;
  
  // brightness data
  int N_bright;
  vector[N_bright] brightness_obs;
  array[N_bright] int<lower=1,upper=N_spp> sp_index_bright;
  
}
transformed data{
  vector[N_bright] log_brightness_obs = log(brightness_obs);
  vector[N_spp] log_island_area = log(island_area);
}
parameters {
  // microhabitat
  real mu_arboreal;
  real<lower=0> sd_arboreal;
  vector[N_spp] arboreal_prob_logit;
  // veg zone
  real mu_arid;
  real<lower=0> sd_arid;
  vector[N_spp] arid_prob_logit;
  // brightness
  real mu_ln_bright;
  real<lower=0> sd_ln_bright_sp;
  real<lower=0> sd_ln_bright_island;
  vector[N_spp] sp_effect_bright;
  vector[N_island] island_effect_bright;
  real slope_arbor;
  real slope_arid;
  real slope_age;
  real slope_area;
  real<lower=0> sigma_bright;
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
  
  // priors on brightness params
  mu_ln_bright ~ normal(8, 2);
  sd_ln_bright_sp ~ exponential(.5);
  sd_ln_bright_island ~ exponential(.5);
  slope_arbor ~ normal(0, .5);
  slope_arid  ~ normal(0, .5);
  slope_age   ~ normal(0, .5);
  slope_area  ~ normal(0, .5);
  sigma_bright ~ exponential(.5);
  // Hierarchical priors for species & island effects on brightness
  sp_effect_bright ~ normal(0, sd_ln_bright_sp);
  island_effect_bright ~ normal(0, sd_ln_bright_island);
  // true brightness -- a SPECIES-level trait
  vector[N_spp] true_ln_bright;
  
  true_ln_bright = mu_ln_bright + 
  sp_effect_bright + 
  island_effect_bright[island_index_spp] + 
  slope_arbor * arboreal_prob + 
  slope_arid * arid_prob + 
  slope_age * island_age +
  slope_area * log_island_area;
  // observed brightness
  log_brightness_obs ~ normal(true_ln_bright[sp_index_bright], sigma_bright);
    
}
generated quantities {
    // true brightness -- a SPECIES-level trait
  vector[N_spp] true_ln_bright;
  
  true_ln_bright = mu_ln_bright + 
  sp_effect_bright + 
  island_effect_bright[island_index_spp] + 
  slope_arbor * arboreal_prob + 
  slope_arid * arid_prob + 
  slope_age * island_age +
  slope_area * log_island_area;
    
}
