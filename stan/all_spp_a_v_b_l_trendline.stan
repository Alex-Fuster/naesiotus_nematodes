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
  
  // load data
  int N_load;
  array[N_load] int load_obs;
  array[N_load] int<lower=1,upper=N_spp> sp_index_load;
  
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
  real slope_arbor_bright;
  real slope_arid_bright;
  real slope_age_bright;
  real slope_area_bright;
  real<lower=0> sigma_bright;
  // parasite load -- no sigma
  real mu_ln_load;
  real<lower=0> sd_ln_load_sp;
  real<lower=0> sd_ln_load_island;
  vector[N_spp] sp_effect_load;
  vector[N_island] island_effect_load;
  real slope_arbor_load;
  real slope_arid_load;
  real slope_age_load;
  real slope_area_load;
  real slope_bright_load;
}
transformed parameters {
  vector[N_spp] arboreal_prob = inv_logit(arboreal_prob_logit);
  vector[N_spp] arid_prob = inv_logit(arid_prob_logit);
}
generated quantities {
  // true brightness -- a SPECIES-level trait
  vector[N_spp] true_ln_bright;
  
  true_ln_bright = mu_ln_bright + 
  sp_effect_bright + 
  island_effect_bright[island_index_spp] + 
  slope_arbor_bright * arboreal_prob + 
  slope_arid_bright * arid_prob + 
  slope_age_bright * island_age +
  slope_area_bright * log_island_area;
  
  vector[N_spp] true_ln_load;
  
  true_ln_load = mu_ln_load + 
  sp_effect_load + 
  island_effect_load[island_index_spp] + 
  slope_arbor_load * arboreal_prob + 
  slope_arid_load   * arid_prob + 
  slope_age_load * island_age +
  slope_area_load * log_island_area + 
  slope_bright_load * true_ln_bright; // hypothesis test!
  // observed load
}