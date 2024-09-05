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

  // Nematode div data
  int N_div;
  array[N_div] int nematode_count;
  array[N_div] int<lower=1, upper=N_spp> sp_index_div;
  array[N_div] int<lower=1, upper=N_island> island_index_div;
  

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

  // Nematode div model parameters
  real mu_div;
  real<lower=0> sd_div_sp;
  real<lower=0> sd_div_island;
  vector[N_spp] sp_effect_div;
  vector[N_island] island_effect_div;
  real slope_bright_div;  // Effect of brightness on nematode div
  real slope_arbor_div;
  real slope_arid_div;
  real slope_age_div;
  real slope_area_div;
  
  real<lower=0> phi_div;  // Dispersion parameter for the negative binomial model
  
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
  slope_arbor_bright ~ normal(0.64, 0.5);  // mostly positive relationship (Kraemer et al. 2019)
  slope_arid_bright ~ normal(0, 0.5); 
  slope_age_bright ~ normal(0, 0.5);
  slope_area_bright ~ normal(0, 0.5);
  sigma_bright ~ exponential(0.5);
  
  // Priors for nematode div model
  mu_div ~ normal(0, 1);
  sd_div_sp ~ exponential(0.5);
  sd_div_island ~ exponential(0.5);
  slope_bright_div ~ normal(0, 0.5);
  slope_arbor_div  ~ normal(0, 0.5); 
  slope_arid_div  ~ normal(0, 0.5); 
  slope_age_div ~ normal(0.64, 0.5); // positive relationship
  slope_area_div ~ normal(0.64, 0.5); // positive relationship
  
  sd_div_sp ~ exponential(0.5);
  sp_effect_div ~ normal(0, sd_div_sp);
  
  sd_div_island ~ exponential(0.5);
  island_effect_div ~ normal(0, sd_div_island);
  
  phi_div ~ exponential(1);


  
  // Likelihood for brightness model
  vector[N_spp] true_ln_bright = mu_ln_bright +
  sp_effect_bright +
  island_effect_bright[island_index_spp] +
  slope_arbor_bright * arboreal_prob +
  slope_arid_bright * arid_prob +
  slope_age_bright * island_age_true +
  slope_area_bright * log_island_area;
  
  log_brightness_obs ~ normal(true_ln_bright[sp_index_bright], sigma_bright);
  
  // Likelihood for nematode div model
  
  vector[N_spp] mean_div = mu_div +
  sp_effect_div +
  island_effect_div[island_index_spp] +
  slope_bright_div * true_ln_bright +  // Use the modeled brightness directly
  slope_arbor_div * arboreal_prob +
  slope_arid_div * arid_prob +
  slope_age_div * island_age_true +
  slope_area_div * log_island_area;
  
  nematode_count ~ neg_binomial_2_log(mean_div[sp_index_div], phi_div);
  
  
}

generated quantities {
  vector[N_spp] true_ln_bright = mu_ln_bright +
                                 sp_effect_bright +
                                 island_effect_bright[island_index_spp] +
                                 slope_arbor_bright * arboreal_prob +
                                 slope_arid_bright * arid_prob +
                                 slope_age_bright * island_age_true +
                                 slope_area_bright * log_island_area;
                                 
  vector[N_div] predicted_nematode_div;
  
  // Generate predicted nematode div using the negative binomial distribution
  for (n in 1:N_div) {
    predicted_nematode_div[n] = neg_binomial_2_log_rng(
      mu_div +
      sp_effect_div[sp_index_div[n]] +
      island_effect_div[island_index_div[n]] +
      slope_bright_div * true_ln_bright[sp_index_div[n]] +  // Use the modeled brightness directly
      slope_arbor_div * arboreal_prob[sp_index_div[n]] +
      slope_arid_div * arid_prob[sp_index_div[n]] +
      slope_age_div * island_age_true[island_index_div[n]] +
      slope_area_div * log_island_area[island_index_div[n]],
      phi_div
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
