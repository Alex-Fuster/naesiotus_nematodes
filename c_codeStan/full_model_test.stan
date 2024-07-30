data {
  int<lower=0> N_nem;               // Number of observations for nematode load
  int<lower=0> N_bright;            // Number of observations for brightness
  int<lower=0> N_habitat;           // Number of observations for habitat
  int<lower=0> N_veg;               // Number of observations for vegetation zone
  
  int<lower=0> N_sp_nem;            // Number of species in the nematode data
  
  array[N_nem] int<lower=0, upper=100> load;               // Nematode load counts
  
  // Load data
  array[N_nem] int<lower=1, upper=N_sp_nem> species_i;     // Species id for load dataset
  
  // Brightness data
  vector[N_bright] brightness;                             // Brightness scores
  array[N_bright] int<lower=1, upper=N_sp_nem> species_k;   // Species id for brightness dataset
  
  // Habitat data
  array[N_habitat] int<lower=0, upper=1> habitat;           // Habitat: 0 for terrestrial, 1 for arboreal
  array[N_habitat] int<lower=1, upper=N_sp_nem> species_h; // Species id for habitat dataset
  
  // Vegetation zone data
  array[N_veg] int<lower=0, upper=1> vegetation;           // Vegetation zone: 0 for humid, 1 for arid
  array[N_veg] int<lower=1, upper=N_sp_nem> species_v;     // Species id for vegetation dataset
  
  // Island age and area data
  vector[N_nem] island_age;            // Island age for each observation
  vector[N_nem] island_area;           // Island area for each observation
}

transformed data {
  vector[N_bright] log_brightness;
  vector[N_nem] log_island_area;
  
  log_brightness = log(brightness);
  log_island_area = log(island_area);  // Log-transform island area
}

parameters {
  real intercept;
  real<upper=0> slope_bright;
  real<upper=0> slope_habitat;
  real<upper=0> slope_vegetation;
  real slope_age;
  real slope_area;
  
  vector[N_sp_nem] mu_bright;          // Mean brightness for each species
  vector<lower=0, upper=1>[N_sp_nem] habitat_prob;    // Habitat probability for each species
  vector<lower=0, upper=1>[N_sp_nem] vegetation_prob; // Vegetation probability for each species
  
  real<lower=0> sigma_bright; // Standard deviation of brightness for each species
}

model {
  vector[N_nem] lambda;
  
  // Priors
  intercept ~ normal(0, 1);
  slope_bright ~ normal(0, 1);
  slope_habitat ~ normal(0, 1);
  slope_vegetation ~ normal(0, 1);
  slope_age ~ normal(0, 1);
  slope_area ~ normal(0, 1);
  
  mu_bright ~ normal(5, 2);
  habitat_prob ~ beta(2, 2);
  vegetation_prob ~ beta(2, 2);
  
  sigma_bright ~ exponential(1);
  
  // Model brightness
  log_brightness ~ normal(mu_bright[species_k], sigma_bright);
  
  // Model habitat
  habitat ~ bernoulli(habitat_prob[species_h]);
  
  // Model vegetation zone
  vegetation ~ bernoulli(vegetation_prob[species_v]);
  
  // Linear predictor combining species-level and island-level predictors
  lambda = intercept + 
           slope_bright * mu_bright[species_k] + 
           slope_habitat * habitat_prob[species_h] + 
           slope_vegetation * vegetation_prob[species_v] + 
           slope_age * island_age + 
           slope_area * log_island_area;
  
  // Model nematode load
  load ~ poisson_log(lambda);
}

generated quantities {
  array[N_bright] real predicted_brightness;  // Predicted brightness
  array[N_habitat] real predicted_habitat;   // Predicted habitat probabilities
  array[N_veg] real predicted_vegetation;    // Predicted vegetation probabilities
  array[N_nem] int predicted_load;           // Predicted nematode load
  
  for (n in 1:N_nem) {
    predicted_load[n] = poisson_log_rng(intercept + 
                                         slope_bright * mu_bright[species_i[n]] + 
                                         slope_habitat * habitat_prob[species_i[n]] + 
                                         slope_vegetation * vegetation_prob[species_i[n]] + 
                                         slope_age * island_age[n] + 
                                         slope_area * log_island_area[n]);
  }
  
  for (b in 1:N_bright) {
    predicted_brightness[b] = exp(normal_rng(mu_bright[species_k[b]], sigma_bright));
  }
  
  for (h in 1:N_habitat) {
    predicted_habitat[h] = bernoulli_rng(habitat_prob[species_h[h]]);
  }
  
  for (v in 1:N_veg) {
    predicted_vegetation[v] = bernoulli_rng(vegetation_prob[species_v[v]]);
  }
}
