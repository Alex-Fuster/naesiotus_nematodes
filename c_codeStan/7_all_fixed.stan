data {
  int<lower = 1> N; // Number of observation
  int<lower = 1> N_Y; // Number of years
  array[N] int<lower = 0, upper = 1> y; // Arrey of our observation (response)
  array[N] int<lower = 1, upper = N_Y> year; // Indices of years
  vector[N] effort; // Vector of effort (explanatory variable)
  vector[N] jj_date; // Vector of observation data (explanatory variable)
  int<lower=0> n_new;
  vector[n_new*N_Y] newdate;
  array[n_new*N_Y] int<lower=1, upper = N_Y> newyear;
}

// The parameters accepted by the model. Our model
// accepts five parameters.
parameters {
  real<lower=0, upper=1> d;
  
  vector[N_Y] log_a; // Value of a
  vector[N_Y] f; // Link with slope of a
  vector<lower=1, upper=365>[N_Y] b1;
  vector<lower=1, upper=365>[N_Y] b2;
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.

transformed parameters{
  vector[N_Y] f_prob = inv_logit(f);
  vector[N_Y] a_abs = exp(log_a);
}

model {
  d ~ beta(2, 5);
  log_a ~  normal(-1, 0.8);
  f ~  normal(1, 0.5);
  b1 ~ normal(150, 7);
  b2 ~ normal(210, 7);
  
  
  y ~ bernoulli((1 - (1 - d)^effort) .* 
                inv_logit(-a_abs[year] .* f_prob[year] .* (jj_date - b1[year])) .*
                inv_logit(a_abs[year] .* (1 - f_prob[year]) .* (jj_date - b2[year])));
}
generated quantities {
  vector[n_new*N_Y] mu_line = inv_logit(-a_abs[newyear] .* f_prob[newyear] .* (newdate - b1[newyear])) .*
                              inv_logit(a_abs[newyear] .* (1 - f_prob[newyear]) .* (newdate - b2[newyear]));
}
