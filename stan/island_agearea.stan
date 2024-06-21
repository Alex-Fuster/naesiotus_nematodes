data {
  int<lower=0> N_nem;                      // Number of observations
  array[N_nem] int<lower=0, upper=100> load;  // Nematode load counts
  vector[N_nem] island_age;               // Island age for each observation
  vector[N_nem] island_area;            // Island area for each observation
}

transformed data {
  vector[N_nem] log_island_area;
  log_island_area = log(island_area);  // Log-transform island area
}

parameters {
  real intercept;
  real slope_age;                      // Slope for island age
  real slope_area;                     // Slope for log-transformed island area
}

model {
  vector[N_nem] lambda;
  
  // Priors
  intercept ~ normal(0, 1);
  slope_age ~ normal(0, 1);
  slope_area ~ normal(0, 1);
  
  // Linear predictor
  lambda = intercept + slope_age * island_age + slope_area * log_island_area;
  
  // Model
  load ~ poisson_log(lambda); 
}

generated quantities {
  array[N_nem] int predicted_load;  // Predicted nematode load
  
  for (n in 1:N_nem) {
    predicted_load[n] = poisson_log_rng(intercept + slope_age * island_age[n] + slope_area * log_island_area[n]);
  }
}
