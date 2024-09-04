data {
  int<lower=0> N_load;                      // Number of observations load dataset
  int<lower=0> N_bright;                      // Number of observations brightness dataset
  int<lower=0> N_spp;                  // Number of unique species
  int<lower=0> N_islands;              // Number of unique islands
  int<lower=0> N_habitat;              // Number of habitat observations
  int<lower=0> N_veg;                  // Number of vegetation zone observations
  
  

  
  array[N_bright] int<lower=1,upper=N_spp> spp_id_bright;   // Species ID for brightness data
  array[N_load] int<lower=1,upper=N_islands> island_id_load;    // Island ID for each observation in the load dataset
  array[N_bright] int<lower=1,upper=N_islands> island_id_bright;    // Island ID for each observation in the brightness dataset
  array[N_habitat] int<lower=1,upper=N_spp> spp_id_hab; // Species ID for habitat data
  array[N_veg] int<lower=1,upper=N_spp> spp_id_veg;     // Species ID for vegetation data
  array[N_load] int<lower=1,upper=N_spp> spp_id_load;     // Species ID for load data
  
  array[N_load] int<lower=0, upper=100> load;               // Nematode load counts
  
  vector[N_bright] brightness;                // Brightness data
  
  array[N_habitat] int habitat_arboreal;   // Total counts of arid habitat for each species
  array[N_habitat] int total_hab;  // Total counts of humid habitat for each species
  
  array[N_veg] int vegetation_arid;    // Total counts of arid vegetation for each species
  array[N_veg] int total_veg;   // Total counts of humid vegetation for each species
  
  vector[N_load] island_age_load;                // Island age for each observation
  vector[N_bright] island_age_bright;                // Island age for each observation
  vector[N_load] island_area_load;               // Island area for each observation
  vector[N_bright] island_area_bright;               // Island area for each observation
  
}

transformed data {
  vector[N_bright] log_brightness = log(brightness);
  vector[N_bright] log_island_area_bright = log(island_area_bright);
  vector[N_load] log_island_area_load = log(island_area_load);

}

parameters {
  
  real intercept;
  real slope_bright;
  real<lower=0> slope_arboreal_bright;         // Slope for arboreal habitat effect on brightness
  real<upper=0> slope_arboreal_load;         // Slope for arboreal habitat effect on load
  real<lower=0> slope_arid_bright;             // Slope for arid vegetation effect on brightness
  real<upper=0> slope_arid_load;             // Slope for arid vegetation effect on load
  real slope_age;                       // Slope for island age effect
  real slope_area;                      // Slope for log island area effect
  vector[N_spp] spp_effects;            // Species random effects
  
  vector[N_islands] z_island_effects;     // Island random effects
  real mu;                              // Overall mean
  
  real<lower=0> sigma_bright;           // Standard deviation for brightness
  real<lower=0> sigma_spp;
  real<lower=0> sigma_islands;
  
  vector[N_spp] arboreal_prob_logit;    // Arboreal habitat probability for each species
  vector[N_spp] arid_prob_logit;        // Arid vegetation probability for each species
  
 
  // Hyperparameters for the species-specific probabilities
  real mu_arboreal;
  real<lower=0> sd_arboreal;
  real mu_arid;
  real<lower=0> sd_arid;
  
}
transformed parameters{
  vector[N_spp] arboreal_prob = inv_logit(arboreal_prob_logit);
  vector[N_spp] arid_prob = inv_logit(arid_prob_logit);
  vector[N_islands] island_effects;
  island_effects = z_island_effects * sigma_islands;
  
}

generated quantities{
  
  vector[N_bright] log_brightness_pred;
  vector[N_bright] brightness_pred;
  vector[N_spp] habitat_arboreal_pred;
  vector[N_spp] vegetation_arid_pred;
  array[N_load] int load_pred;
  
  // vector[N_islands] island_effects;
  // island_effects = z_island_effects * sigma_islands;
  
  for (i in 1:N_bright) {
    log_brightness_pred[i] = normal_rng(
      mu +
      spp_effects[spp_id_bright[i]] +
      island_effects[island_id_bright[i]] +
      slope_arboreal_bright * arboreal_prob[spp_id_bright[i]] +
      slope_arid_bright * arid_prob[spp_id_bright[i]] +
      slope_age * island_age_bright[island_id_bright[i]] +
      slope_area * log_island_area_bright[island_id_bright[i]],
      sigma_bright);
  }
  
  brightness_pred = exp(log_brightness_pred);
  
  for (s in 1:N_habitat) {
    habitat_arboreal_pred[s] = binomial_rng(total_hab[s], arboreal_prob[s]);
  }
  
  for (s in 1:N_veg) {
    vegetation_arid_pred[s] = binomial_rng(total_veg[s], arid_prob[s]);
  }
  
  
  for (i in 1:N_load) {
    load_pred[i] = poisson_log_rng(
      intercept + 
    slope_bright * log_brightness[spp_id_load[i]] + 
    slope_arboreal_load * arboreal_prob[spp_id_load[i]] +
    slope_arid_load * arid_prob[spp_id_load[i]] +
    slope_age * island_age_load[island_id_load[i]] +
    slope_area * log_island_area_load[island_id_load[i]]);

  }
  

}
