data {
  int N_spp;  // Number of species (47)
  int N_island;  // Number of islands
  array[N_spp] int<lower=1, upper=N_spp> sp_id;
  vector[N_spp] island_area;
 // vector[N_spp] island_age;
   vector[N_spp] min_emergence; // Minimum estimated age
  vector[N_spp] max_emergence; // Maximum estimated age
  array[N_spp] int<lower=1,upper=N_island> island_index_spp;
  
  // Habitat data
  int N_habitat;
  array[N_habitat] int habitat_arboreal;
  array[N_habitat] int total_hab;
  array[N_habitat] int<lower=1,upper=N_spp> sp_index_hab;
  
  // Vegetation zone data
  int N_veg;
  array[N_veg] int veg_arid;
  array[N_veg] int total_veg;
  array[N_veg] int<lower=1,upper=N_spp> sp_index_veg;
  
  // Brightness data
  int N_bright;
  vector[N_bright] brightness_obs;
  array[N_bright] int<lower=1,upper=N_spp> sp_index_bright;

  // Nematode load data
  int N_load;
  array[N_load] int nematode_count;
  array[N_load] int<lower=1, upper=N_spp> sp_index_load;
  array[N_load] int<lower=1, upper=N_island> island_index_load;
  

}

transformed data {
  vector[N_bright] log_brightness_obs = log(brightness_obs);
  vector[N_spp] log_island_area = log(island_area);
}

parameters {
  
  // Latent parameter for island age
  vector<lower=min_emergence, upper=max_emergence>[N_spp] island_age_true;  // True island ages
  
  // microhabitat
  real mu_arboreal;
  real<lower=0> sd_arboreal;
  vector[N_spp] arboreal_prob_logit;
  // veg zone
  real mu_arid;
  real<lower=0> sd_arid;
  vector[N_spp] arid_prob_logit;
  
  // Brightness model parameters
  real mu_ln_bright;
  real<lower=0> sd_ln_bright_sp;
  real<lower=0> sd_ln_bright_island;
  vector[N_spp] sp_effect_bright;
  vector[N_island] island_effect_bright;
  real slope_arbor_bright;
  real slope_arid_bright;
  real slope_age_bright;
  real slope_area_bright;
  real<lower=0> sigma_bright;

  // Nematode load model parameters
  real mu_load;
  real<lower=0> sd_load_sp;
  real<lower=0> sd_load_island;
  vector[N_spp] sp_effect_load;
  vector[N_island] island_effect_load;
  real slope_bright_load;  // Effect of brightness on nematode load
  real slope_arbor_load;
  real slope_arid_load;
  real slope_age_load;
  real slope_area_load;
  
  real<lower=0> phi_load;  // Dispersion parameter for the negative binomial model
  
}

transformed parameters {
  
  vector[N_spp] arboreal_prob = inv_logit(arboreal_prob_logit);
  vector[N_spp] arid_prob = inv_logit(arid_prob_logit);
}

model {
  
   // Priors for latent island age
  island_age_true ~ uniform(min_emergence, max_emergence); // Uniform prior between min and max age
  
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
  
  
  
  // Priors for brightness model
  mu_ln_bright ~ normal(8, 2);
  sd_ln_bright_sp ~ exponential(0.5);
  sd_ln_bright_island ~ exponential(0.5);
  sp_effect_bright ~ normal(0, sd_ln_bright_sp);
  island_effect_bright ~ normal(0, sd_ln_bright_island);
  slope_arbor_bright ~ normal(0, 0.5); 
  slope_arid_bright ~ normal(0, 0.5);
  slope_age_bright ~ normal(0, 0.5); 
  slope_area_bright ~ normal(0, 0.5);
  sigma_bright ~ exponential(0.5);
  
  // Priors for nematode load model
  mu_load ~ normal(0, 1);

  slope_bright_load ~ normal(0, 0.5);
  slope_arbor_load  ~ normal(0, 0.5); 
  slope_arid_load  ~ normal(0, 0.5); 
  slope_age_load ~ normal(0, 0.5);
  slope_area_load ~ normal(0, 0.5);
  
  sd_load_sp ~ exponential(0.5);
  sp_effect_load ~ normal(0, sd_load_sp);
  
  sd_load_island ~ exponential(0.5);
  island_effect_load ~ normal(0, sd_load_island);
  
  phi_load ~ exponential(1);


  
  // Likelihood for brightness model
  vector[N_spp] true_ln_bright = mu_ln_bright +
  sp_effect_bright +
  island_effect_bright[island_index_spp] +
  slope_arbor_bright * arboreal_prob +
  slope_arid_bright * arid_prob +
  slope_age_bright * island_age_true +
  slope_area_bright * log_island_area;
  
  log_brightness_obs ~ normal(true_ln_bright[sp_index_bright], sigma_bright);
  
  // Likelihood for nematode load model
  
  vector[N_spp] mean_load = mu_load +
  sp_effect_load +
  island_effect_load[island_index_spp] +
  slope_bright_load * true_ln_bright +  // Use the modeled brightness directly
  slope_arbor_load * arboreal_prob +
  slope_arid_load * arid_prob +
  slope_age_load * island_age_true +
  slope_area_load * log_island_area;
  
  nematode_count ~ neg_binomial_2_log(mean_load[sp_index_load], phi_load);
  
  
}

generated quantities {
  vector[N_spp] true_ln_bright = mu_ln_bright +
                                 sp_effect_bright +
                                 island_effect_bright[island_index_spp] +
                                 slope_arbor_bright * arboreal_prob +
                                 slope_arid_bright * arid_prob +
                                 slope_age_bright * island_age_true +
                                 slope_area_bright * log_island_area;
                                 
  vector[N_load] predicted_nematode_load;
  
  // Generate predicted nematode load using the negative binomial distribution
  for (n in 1:N_load) {
    predicted_nematode_load[n] = neg_binomial_2_log_rng(
      mu_load +
      sp_effect_load[sp_index_load[n]] +
      island_effect_load[island_index_load[n]] +
      slope_bright_load * true_ln_bright[sp_index_load[n]] +  // Use the modeled brightness directly
      slope_arbor_load * arboreal_prob[sp_index_load[n]] +
      slope_arid_load * arid_prob[sp_index_load[n]] +
      slope_age_load * island_age_true[island_index_load[n]] +
      slope_area_load * log_island_area[island_index_load[n]],
      phi_load
    );
  }
  
  // Generated quantities for habitat and vegetation zone
  vector[N_habitat] predicted_habitat_arboreal;  // Predicted number of species in arboreal habitat
  vector[N_veg] predicted_veg_arid;  // Predicted number of species in arid vegetation zone
  
  // Generate predicted counts for habitat (arboreal)
  for (h in 1:N_habitat) {
    predicted_habitat_arboreal[h] = binomial_rng(total_hab[h], inv_logit(arboreal_prob_logit[sp_index_hab[h]]));
  }
  
  // Generate predicted counts for vegetation zone (arid)
  for (v in 1:N_veg) {
    predicted_veg_arid[v] = binomial_rng(total_veg[v], inv_logit(arid_prob_logit[sp_index_veg[v]]));
  }
  
  
}
