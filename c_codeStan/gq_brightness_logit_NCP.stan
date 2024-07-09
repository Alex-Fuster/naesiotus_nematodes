data {
  int<lower=0> N;                      // Number of observations
  int<lower=0> N_spp;                  // Number of unique species
  int<lower=0> N_islands;              // Number of unique islands
  int<lower=0> N_habitat;              // Number of habitat observations
  int<lower=0> N_veg;                  // Number of vegetation zone observations
  
  array[N] int<lower=1,upper=N_spp> spp_id_bright;   // Species ID for brightness data
  array[N] int<lower=1,upper=N_islands> island_id;    // Island ID for each observation
  array[N_habitat] int<lower=1,upper=N_spp> spp_id_hab; // Species ID for habitat data
  array[N_veg] int<lower=1,upper=N_spp> spp_id_veg;     // Species ID for vegetation data
  
  vector[N] brightness;                // Brightness data
  vector[N] island_age;                // Island age for each observation
  vector[N] island_area;               // Island area for each observation
  
  array[N_habitat] int habitat_arboreal;   // Total counts of arid habitat for each species
  array[N_habitat] int total_hab;  // Total counts of humid habitat for each species
  
  array[N_veg] int vegetation_arid;    // Total counts of arid vegetation for each species
  array[N_veg] int total_veg;   // Total counts of humid vegetation for each species
}

transformed data {
  vector[N] log_brightness = log(brightness);
  vector[N] log_island_area = log(island_area);

}

parameters {

  real slope_arboreal;         // Slope for arboreal habitat effect
  real slope_arid;             // Slope for arid vegetation effect
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
generated quantities {
  
  vector[N] log_brightness_pred;
  vector[N] brightness_pred;
  vector[N_spp] habitat_arboreal_pred;
  vector[N_spp] vegetation_arid_pred;
  
 // vector[N_islands] island_effects;
 // island_effects = z_island_effects * sigma_islands;
  
  for (i in 1:N) {
    log_brightness_pred[i] = normal_rng(
      mu +
      spp_effects[spp_id_bright[i]] +
      island_effects[island_id[i]] +
      slope_arboreal * arboreal_prob[spp_id_bright[i]] +
      slope_arid * arid_prob[spp_id_bright[i]] +
      slope_age * island_age[island_id[i]] +
      slope_area * log_island_area[island_id[i]],
      sigma_bright);
  }
  
  brightness_pred = exp(log_brightness_pred);
  
  for (s in 1:N_spp) {
    habitat_arboreal_pred[s] = binomial_rng(total_hab[s], arboreal_prob[s]);
    vegetation_arid_pred[s] = binomial_rng(total_veg[s], arid_prob[s]);
  }
}