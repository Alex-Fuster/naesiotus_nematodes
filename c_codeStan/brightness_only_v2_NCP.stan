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
  
  vector<lower=0, upper=1>[N_spp] arboreal_prob;    // Arboreal habitat probability for each species
  vector<lower=0, upper=1>[N_spp] arid_prob;        // Arid vegetation probability for each species
  
 
  // Hyperparameters for the species-specific probabilities
  real<lower=0> alpha_arboreal;
  real<lower=0> beta_arboreal;
  real<lower=0> alpha_arid;
  real<lower=0> beta_arid;
  
}

model {
  // Priors
  slope_arboreal ~ normal(0, 1);
  slope_arid ~ normal(0, 1);
  slope_age ~ normal(0, 1);
  slope_area ~ normal(0, 1);
  mu ~ normal(5, 2);
  arboreal_prob ~ beta(2, 2);
  arid_prob ~ beta(2, 2);
  sigma_bright ~ exponential(1);
  
  // Species and island effects
  spp_effects ~ normal(0, sigma_spp);
  
  // non-centered
  vector[N_islands] island_effects;
  z_island_effects ~ std_normal();
  island_effects = z_island_effects * sigma_islands;
  
  
  // Hierarchical priors for the probabilities
  arboreal_prob ~ beta(alpha_arboreal, beta_arboreal);
  arid_prob ~ beta(alpha_arid, beta_arid);
  
  
  // Likelihood
  
  habitat_arboreal ~ binomial(total_hab, arboreal_prob[spp_id_hab]);
  vegetation_arid ~ binomial(total_veg, arid_prob[spp_id_veg]);
  
  log_brightness ~ normal(
    mu +
    spp_effects[spp_id_bright] +
    island_effects[island_id] +
    slope_arboreal * arboreal_prob[spp_id_bright] +
    slope_arid * arid_prob[spp_id_bright] +
    slope_age * island_age[island_id] +
    slope_area * log_island_area[island_id],
    sigma_bright);
    
    
}
