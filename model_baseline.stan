data {
    // Number of total observations
    int<lower=1> N;
    
    // Number of unique enhancers (excluding control)
    int<lower=1> E;
    
    // Number of unique guides (excluding control)
    int<lower=1> G;
    
    // Observed data
    real expression[N];            // Expression values
    int<lower=-1> enhancer_idx[N]; // Enhancer index for each observation (-1 for control)
    int<lower=-1> guide_idx[N];    // Guide index for each observation (-1 for control)
}

transformed data {
    // Calculate control mean and standard deviation
    int n_control = 0;
    real control_sum = 0;
    real control_mean;
    
    for (n in 1:N) {
        if (enhancer_idx[n] == 0) {
            n_control += 1;
            control_sum += expression[n];
        }
    }
    control_mean = control_sum / n_control;

}

parameters {    
    // Observation noise standard deviation
    real<lower=0> sigma;
    
    // Enhancer effect sizes (between 0 and 1)
    vector<lower=0, upper=1>[E] b;
    
    // Global probability of guide being functional
    real<lower=0, upper=1> p_functional;
}

model {
    // Priors
    sigma ~ cauchy(0, 2.5);                     // Weakly informative prior for sigma
    b ~ beta(1, 1);                             // Uniform prior on effect sizes

    // Likelihood
    for (n in 1:N) {
        if (enhancer_idx[n] == 0) {
            // Control case
            expression[n] ~ normal(control_mean, sigma);
        } else {
            // Treatment case
            expression[n] ~ normal(control_mean * b[enhancer_idx[n]], sigma);
        }
    }
}
