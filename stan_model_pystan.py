import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
import stan

# Step 1: Load the data
file_path = "data1.txt"  # Update with your actual file path
data = pd.read_csv(file_path, sep="\t")

# Step 2: Divide data into separate files based on geneID
input_folder = "input"  # Folder to save the divided files

# Create the output folder if it doesn't exist
if not os.path.exists(input_folder):
    os.makedirs(input_folder)

# Group data by geneID and prepare data for Stan model
stan_data_list = []

for gene_id, group in data.groupby('geneID'):
    stan_data = {
        'N': len(group),
        'E': group['enhancerID'].nunique() - 1,  # Excluding control (-1)
        'G': group['guideID'].nunique() - 1,     # Excluding control (-1)
        'expression': group['expression'].astype(float).tolist(),
        'enhancer_idx': group['enhancerID'].tolist(),
        'guide_idx': group['guideID'].tolist()
    }
    stan_data_list.append(stan_data)

stan_model_code = """
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
        if (enhancer_idx[n] == -1) {
            n_control += 1;
            control_sum += expression[n];
        }
    }
    
    control_mean = control_sum / n_control;

    // Calculate control standard deviation
    real control_ssq = 0;
    real control_sigma;
    
    for (n in 1:N) {
        if (enhancer_idx[n] == -1) {
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
    b ~ beta(1, 1);                            // Uniform prior on effect sizes
    p_functional ~ beta(2, 2);                  // Prior for global guide functionality probability

    // Likelihood
    for (n in 1:N) {
        if (enhancer_idx[n] == -1) {
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
"""


output_folder = "output"  # Folder to save the summary files

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

num_gene = 1
# Fit the model to your data
for i in range(num_gene):
    # print(stan_data_list[i])
    # Compile the model
    stan_model = stan.build(stan_model_code, data=stan_data_list[i])
    fit = stan_model.sample(num_chains=4, num_samples=1000, data=stan_data_list[i])

    # save the summary
    summary = fit.summary()
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(output_folder, f"summary_gene_{i}.csv"))
