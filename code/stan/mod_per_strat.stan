// Model description:
// Perform Bayesian hierarchical Nowcast in S different strata
// per stratum s: model log(lambda_{t,s}) = N(log(lambda_{t-1, s}), sigma_s), 
// with strata-specific sd
// model observed case counts n_{t,d,s} ~ NB(lambda[t, s]*p_{t,d}, phi_{s}), 
// with strata-specific over-dispersion
// delay distribution: Discrete time-hazard model with week-day effects + 
// strata-specific dummy per week
data {
  // data
  int T;              // Number of rows in reporting triangle 
  int D;              // Maximum delay and number of columns -1 of reporting triangle'
  int r[T, D+1];   // Reporting triangle (Including zero)
  int k_wd_haz; // number weekday-effects discrete-time hazard model
  matrix[T, k_wd_haz] W_wd[D];  // Design matrix for discrete hazard model
  int n_week;
  int week_ind[T];
  // prior parameter
  vector[D+1] alpha;  // Parameters of Dirichlet prior for baseline delay distribution
}



parameters {
  simplex[D+1] p_bl_pr; // delay probabilities
  vector[T] logLambda_raw;   // Standard normal random-walk changes (scaled by sigma)
  // hospitalization model
  real<lower=0, upper=1> sigma; // Variance parameter for random walk
  // reporting model
  // week-day effect
  vector[k_wd_haz] beta_wd_haz;       
  // week effect
  vector[n_week] beta_week_haz;       
  // Hyperprior
  real<lower=0> sd_beta_wd_haz;
  real<lower=0> sd_beta_week_haz;

  // data model
  real<lower=0, upper=1> reciprocal_phi;   // dispersion parameter: var=mu+reciprocal_phi*mu^2
  real<lower=0, upper=1> reciprocal_phi_d0;   // dispersion parameter d=0: var=mu+reciprocal_phi*mu^2
  real<lower=0, upper=1> reciprocal_phi_d1;   // dispersion parameter d=0: var=mu+reciprocal_phi*mu^2
  real<lower=0, upper=1> reciprocal_phi_d2;   // dispersion parameter d=0: var=mu+reciprocal_phi*mu^
  real<lower=0, upper=1> reciprocal_phi_dmax;   // dispersion parameter d=0: var=mu+reciprocal_phi*mu^2

}

transformed parameters {
  // hospitalization model
  vector[T] logLambda; // expected number of cases at time t in strata s

  // reporting model
  vector[D] gamma;   // Discrete hazard model intercept
  matrix[T, D+1] h;        // Discrete hazard w.r.t time and delay
  matrix[T, D+1] p;        // Reporting probability w.r.t. time and delay
  
  // data model
  vector[D+1] phi;
  // Discrete hazard model
  gamma = logit((p_bl_pr ./ cumulative_sum(p_bl_pr[(D+1):1])[(D+1):1])[1:D]);
  
  for (d in 1:D){
    h[, d] = inv_logit(gamma[d] + W_wd[d]*beta_wd_haz + sd_beta_week_haz*beta_week_haz[week_ind]);
    if (d==1) {
        p[, d] = h[, d];
      } else {      
        p[, d] = h[, d] .* (1 - (p[, 1:(d-1)] * rep_vector(1, d-1)));
      }
  }
  h[, D+1] = rep_vector(1, T);
  p[, D+1] = 1 -  (p[, 1:D] * rep_vector(1, D));
  // log-lambda
  logLambda[1] = logLambda_raw[1];
  for(t in 2:T) {
    logLambda[t] = logLambda[t-1] + sigma*logLambda_raw[t]; // Derive logLambda from non-centered parametrization
  }
  // Overdispersion
  phi[1] = 1 / reciprocal_phi_d0;
  phi[2] = 1 / reciprocal_phi_d1;
  phi[3] = 1 / reciprocal_phi_d2;
  phi[4:D] = rep_vector(1/reciprocal_phi, D-3);
  phi[D+1] =  1 / reciprocal_phi_dmax;
}

model {
  // Priors
  // hospitalization model
  sigma ~ normal(0, .25); // scale of S first-order random walk(s)

  // reporting delay
  // Hyper-prior
  p_bl_pr ~ dirichlet(alpha);
  sd_beta_wd_haz ~ normal(0, 0.5);
  sd_beta_week_haz ~ normal(0, 0.5);
  // Prior
  beta_wd_haz ~ normal(0, sd_beta_wd_haz);
  beta_week_haz ~ std_normal();
  
  // log-Lambda

    logLambda_raw[1] ~ normal(log(max(sum(r[1,]),1)), .1); // logLambda_raw[1] = logLambda[1], prior centered around actual observation
    logLambda_raw[2:T] ~ std_normal(); // standard normal random-walk (non-centered parametrization), scaled by sigma

  // Model for observed counts
  
  for (t in 1:T) {
    r[t, 1:min((T - t)+1, (D+1))] ~ neg_binomial_2(exp(logLambda[t]) * 
    to_vector(p[t, 1:min((T - t)+1, (D+1))]), phi[1:min((T - t)+1, (D+1))]);
  }
}

generated quantities {
  // Define matrix with relevant case counts n_{t,D} 
  // (consider all days with at least one unobserved reporting number (all t with t+D>T))
  int n[D, D+1];
  // Define vector N with sum of event counts over all delays 
  int N[D];
  int N_7d[D];

  // Fill entries of n with observed case counts when available, or sample from 
  // posterior predictive if unaivalable
  for (t in (T-D+1):T){
    for (d in 1:(D+1)){
      if (t + (d-1) <= T) { // Delay d starts with zero delay
      n[(t-T+D), d] = r[t, d];
      } else {
        n[(t-T+D), d] = neg_binomial_2_rng(exp(logLambda[t]) * p[t, d], 
        phi[d]);
      }
    }
    
    // Sum over all columns of n_td to obtain posterior predictive of overall number of hospitalizations per day
    N[t-T+D] = sum(n[t-T+D,]);
    // Sum over last 7 days to obtain posterior predictive of 7day hospitalization
    if (t-T+D>=7) {
      N_7d[t-T+D] = sum(N[(t-T+D-6):(t-T+D)]);
    } else {
      N_7d[t-T+D] = sum(N[1:(t-T+D)]) + sum(to_array_1d(r[(t-6):(T-D),]));
    }
  }
}
