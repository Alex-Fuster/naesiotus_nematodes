data {
  int N_spp;  // Number of species (47)
  int N_island;  // Number of islands
  array[N_spp] int<lower=1, upper=N_spp> sp_id;

  // Island-level data
  vector[N_island] min_emergence; // Minimum estimated age for each island
  vector[N_island] max_emergence; // Maximum estimated age for each island
  vector[N_island] island_area; // Area for each island
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

  // Nematode diversity data
  int N_div;
  vector[N_div] faith_pd;  // Observed Faith PD
  array[N_div] int<lower=1, upper=N_spp> sp_index_div;
  array[N_div] int<lower=1, upper=N_island> island_index_div;
}

transformed data {
  vector[N_bright] log_brightness_obs = log(brightness_obs);
  vector[N_island] log_island_area = log(island_area);
}

parameters {
  // Latent parameter for island age
  vector<lower=min_emergence, upper=max_emergence>[N_island] island_age_true;  // True island ages
  
  // Microhabitat
  real mu_arboreal;
  real<lower=0> sd_arboreal;
  vector[N_spp] arboreal_prob_logit;
  
  // Vegetation zone
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

  // Faith PD model parameters
  real mu_div;
  real<lower=0> sd_div_sp;
  real<lower=0> sd_div_island;
  vector[N_spp] sp_effect_div;
  vector[N_island] island_effect_div;
  real slope_bright_div;  // Effect of brightness on nematode diversity
  real slope_arbor_div;
  real slope_arid_div;
  real slope_age_div;
  real slope_area_div;
  real<lower=0> sigma_div;  // Standard deviation for Faith PD log-normal model
}

transformed parameters {
  vector[N_spp] arboreal_prob = inv_logit(arboreal_prob_logit);
  vector[N_spp] arid_prob = inv_logit(arid_prob_logit);
}

model {
  // Priors for latent island age
  island_age_true ~ uniform(min_emergence, max_emergence); // Uniform prior between min and max age
  
  // Priors for microhabitat
  mu_arboreal ~ std_normal();
  sd_arboreal ~ exponential(1);
  arboreal_prob_logit ~ normal(mu_arboreal, sd_arboreal);
  habitat_arboreal ~ binomial_logit(total_hab, arboreal_prob_logit[sp_index_hab]);
  
  // Priors for vegetation zone
  mu_arid ~ std_normal();
  sd_arid ~ exponential(1);
  arid_prob_logit ~ normal(mu_arid, sd_arid);
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
  
  // Priors for Faith PD model
  mu_div ~ normal(0, 1);
  sd_div_sp ~ exponential(0.5);
  sp_effect_div ~ normal(0, sd_div_sp);
  sd_div_island ~ exponential(0.5);
  island_effect_div ~ normal(0, sd_div_island);
  slope_bright_div ~ normal(0, 0.5);
  slope_arbor_div ~ normal(0, 0.5);
  slope_arid_div ~ normal(0, 0.5);
  slope_age_div ~ normal(0, 0.5);
  slope_area_div ~ normal(0, 0.5);
  sigma_div ~ exponential(0.5);

  // Brightness likelihood
  vector[N_spp] true_ln_bright = mu_ln_bright +
    sp_effect_bright +
    island_effect_bright[island_index_spp] +
    slope_arbor_bright * arboreal_prob +
    slope_arid_bright * arid_prob +
    slope_age_bright * island_age_true[island_index_spp] +
    slope_area_bright * log_island_area[island_index_spp];
  log_brightness_obs ~ normal(true_ln_bright[sp_index_bright], sigma_bright);

  // Faith PD likelihood
  vector[N_div] mean_div = mu_div +
    sp_effect_div[sp_index_div] +
    island_effect_div[island_index_div] +
    slope_bright_div * true_ln_bright[sp_index_div] +
    slope_arbor_div * arboreal_prob[sp_index_div] +
    slope_arid_div * arid_prob[sp_index_div] +
    slope_age_div * island_age_true[island_index_div] +
    slope_area_div * log_island_area[island_index_div];
  faith_pd ~ lognormal(mean_div, sigma_div);
}

generated quantities {
  
  
  vector[N_spp] true_ln_bright = mu_ln_bright +
  sp_effect_bright +
  island_effect_bright[island_index_spp] +
  slope_arbor_bright * arboreal_prob +
  slope_arid_bright * arid_prob +
  slope_age_bright * island_age_true[island_index_spp] +
  slope_area_bright * log_island_area[island_index_spp];
  
  
  
  // Predicted Faith PD
  vector[N_div] predicted_faith_pd;
  for (n in 1:N_div) {
    predicted_faith_pd[n] = lognormal_rng(mu_div +
      sp_effect_div[sp_index_div[n]] +
      island_effect_div[island_index_div[n]] +
      slope_bright_div * true_ln_bright[sp_index_div[n]] +
      slope_arbor_div * arboreal_prob[sp_index_div[n]] +
      slope_arid_div * arid_prob[sp_index_div[n]] +
      slope_age_div * island_age_true[island_index_div[n]] +
      slope_area_div * log_island_area[island_index_div[n]],
      sigma_div);
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

