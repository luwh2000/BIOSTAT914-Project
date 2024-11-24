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
    // Calculate number of control cases and control mean
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

    // Calculate control standard deviation
    real control_ssq = 0;
    real control_sigma;
    
    for (n in 1:N) {
        if (enhancer_idx[n] == 0) {
            control_ssq += square(expression[n] - control_mean);
        }
    }
    
    control_sigma = sqrt(control_ssq / (n_control - 1));
}

parameters {
    // Control expression mean (X)
    real X;
    
    // Observation noise standard deviation
    real<lower=0> sigma;
    
    // Enhancer effect sizes (between 0 and 1)
    vector<lower=0, upper=1>[E] b;
    
    // Global probability of guide being functional
    real<lower=0, upper=1> p_functional;
}

model {
    // Priors
    X ~ normal(control_mean, control_sigma);    // Prior for control mean expression
    sigma ~ cauchy(0, 2.5);                     // Weakly informative prior for sigma
    b ~ beta(1, 0.5);                            // Has a tail at 1
    p_functional ~ beta(2, 2);                  // Prior for global guide functionality probability

    // Likelihood
    for (n in 1:N) {
        if (enhancer_idx[n] == 0) {
            // Control case
            expression[n] ~ normal(X, sigma);
        } else {
            // Treatment case - marginalized over guide functionality
            // Instead of using discrete guide_functional parameter, we directly use p_functional
            target += log_mix(p_functional,
                            normal_lpdf(expression[n] | X * b[enhancer_idx[n]], sigma),
                            normal_lpdf(expression[n] | X, sigma));
        }
    }
}

// generated quantities {
//     // Generate predicted values for model checking
//     vector[N] y_pred;
    
//     for (n in 1:N) {
//         if (enhancer_idx[n] == -1) {
//             y_pred[n] = normal_rng(X, sigma);
//         } else {
//             // For predictions, we can simulate the discrete choice
//             int is_functional = bernoulli_rng(p_functional);
            
//             if (is_functional) {
//                 y_pred[n] = normal_rng(X * b[enhancer_idx[n]], sigma);
//             } else {
//                 y_pred[n] = normal_rng(X, sigma);
//             }
//         }
//     }
// }
