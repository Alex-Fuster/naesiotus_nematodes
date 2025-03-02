data {
  int N_spp;  // Number of species (47)
  int N_island;  // Number of islands
  array[N_spp] int<lower=1, upper=N_spp> sp_id;

  // Island-level data
  vector[N_island] min_emergence;
  vector[N_island] max_emergence;
  vector[N_island] island_area;
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

  // Mismatch data (Replacing brightness)
  int N_mismatch;
  vector[N_mismatch] mismatch_obs;
  array[N_mismatch] int<lower=1,upper=N_spp> sp_index_mismatch;

  // Nematode load data
  int N_load;
  array[N_load] int nematode_count;
  array[N_load] int<lower=1, upper=N_spp> sp_index_load;
  array[N_load] int<lower=1, upper=N_island> island_index_load;
}

transformed data {
  vector[N_island] log_island_area = log(island_area);
}

parameters {
  // Latent parameter for island age
  vector<lower=min_emergence, upper=max_emergence>[N_island] island_age_true;

  // Microhabitat
  real mu_arboreal;
  real<lower=0> sd_arboreal;
  vector[N_spp] arboreal_prob_logit;

  // Vegetation zone
  real mu_arid;
  real<lower=0> sd_arid;
  vector[N_spp] arid_prob_logit;

  // Mismatch model parameters (Skew-Normal)
  real mu_mismatch;
  real<lower=0> sigma_mismatch;
  real alpha_mismatch;
  vector[N_spp] sp_effect_mismatch;
  vector[N_island] island_effect_mismatch;
  real slope_arbor_mismatch;
  real slope_arid_mismatch;
  real slope_age_mismatch;
  real slope_area_mismatch;

  // Nematode load model parameters
  real mu_load;
  real<lower=0> sd_load_sp;
  real<lower=0> sd_load_island;
  vector[N_spp] sp_effect_load;
  vector[N_island] island_effect_load;
  real slope_mismatch_load;
  real slope_arbor_load;
  real slope_arid_load;
  real slope_age_load;
  real slope_area_load;
  real<lower=0> phi_load;
}

transformed parameters {
  vector[N_spp] arboreal_prob = inv_logit(arboreal_prob_logit);
  vector[N_spp] arid_prob = inv_logit(arid_prob_logit);
}

model {
  // Priors for latent island age
  island_age_true ~ uniform(min_emergence, max_emergence);

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

  // Priors for mismatch (Skew-Normal)
  mu_mismatch ~ normal(0, 5);
  sigma_mismatch ~ exponential(0.5);
  alpha_mismatch ~ normal(0, 2);
  sp_effect_mismatch ~ normal(0, sigma_mismatch);
  island_effect_mismatch ~ normal(0, sigma_mismatch);
  slope_arbor_mismatch ~ normal(0, 0.5);
  slope_arid_mismatch ~ normal(0, 0.5);
  slope_age_mismatch ~ normal(0, 0.5);
  slope_area_mismatch ~ normal(0, 0.5);

  // Likelihood for mismatch (Skew-Normal)
  vector[N_spp] true_mismatch = mu_mismatch +
    sp_effect_mismatch +
    island_effect_mismatch[island_index_spp] +
    slope_arbor_mismatch * arboreal_prob +
    slope_arid_mismatch * arid_prob +
    slope_age_mismatch * island_age_true[island_index_spp] +
    slope_area_mismatch * log_island_area[island_index_spp];

  mismatch_obs ~ skew_normal(true_mismatch[sp_index_mismatch], sigma_mismatch, alpha_mismatch);

  // Priors for nematode load model
  mu_load ~ normal(0, 1);
  sd_load_sp ~ exponential(0.5);
  sp_effect_load ~ normal(0, sd_load_sp);
  sd_load_island ~ exponential(0.5);
  island_effect_load ~ normal(0, sd_load_island);
  phi_load ~ exponential(1);
  slope_mismatch_load ~ normal(0, 0.5);
  slope_arbor_load ~ normal(0, 0.5);
  slope_arid_load ~ normal(0, 0.5);
  slope_age_load ~ normal(0, 0.5);
  slope_area_load ~ normal(0, 0.5);

  // Likelihood for nematode load model
  vector[N_spp] mean_load = mu_load +
    sp_effect_load +
    island_effect_load[island_index_spp] +
    slope_mismatch_load * true_mismatch +
    slope_arbor_load * arboreal_prob +
    slope_arid_load * arid_prob +
    slope_age_load * island_age_true[island_index_spp] +
    slope_area_load * log_island_area[island_index_spp];

  nematode_count ~ neg_binomial_2_log(mean_load[sp_index_load], phi_load);
}

generated quantities {
  // Compute posterior predictions for mismatch
  vector[N_spp] true_mismatch = mu_mismatch +
    sp_effect_mismatch +
    island_effect_mismatch[island_index_spp] +
    slope_arbor_mismatch * arboreal_prob +
    slope_arid_mismatch * arid_prob +
    slope_age_mismatch * island_age_true[island_index_spp] +
    slope_area_mismatch * log_island_area[island_index_spp];

  // Compute posterior predictions for nematode load
  vector[N_spp] log_avg_predicted_load = mu_load +
    sp_effect_load +
    island_effect_load[island_index_spp] +
    slope_mismatch_load * true_mismatch +
    slope_arbor_load * arboreal_prob +
    slope_arid_load * arid_prob +
    slope_age_load * island_age_true[island_index_spp] +
    slope_area_load * log_island_area[island_index_spp];

  vector[N_spp] avg_predicted_load = exp(log_avg_predicted_load);

  // Predict individual nematode loads
  vector[N_load] predicted_nematode_load;
  for (n in 1:N_load) {
    predicted_nematode_load[n] = neg_binomial_2_log_rng(
      log_avg_predicted_load[sp_index_load[n]],
      phi_load
    );
  }

  // Predict habitat and vegetation zones
  vector[N_habitat] predicted_habitat_arboreal;
  vector[N_veg] predicted_veg_arid;

  for (h in 1:N_habitat) {
    predicted_habitat_arboreal[h] = binomial_rng(total_hab[h], inv_logit(arboreal_prob_logit[sp_index_hab[h]]));
  }

  for (v in 1:N_veg) {
    predicted_veg_arid[v] = binomial_rng(total_veg[v], inv_logit(arid_prob_logit[sp_index_veg[v]]));
  }
}
