import pandas as pd
import os
import json

# Step 1: Load the data
file_path = "data1.txt"  # Update with your actual file path
data = pd.read_csv(file_path, sep="\t")

# Step 2: Divide data into separate files based on geneID
output_folder = "input"  # Folder to save the divided files

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Group data by geneID and save each group to a separate JSON file
for gene_id, group in data.groupby('geneID'):
    output_file = os.path.join(output_folder, f"gene_{gene_id}.json")
    
    # Map enhancerID and guideID to start from 1
    enhancer_mapping = {old_id: new_id for new_id, old_id in enumerate(sorted(group['enhancerID'].unique()), start=0)}
    guide_mapping = {old_id: new_id for new_id, old_id in enumerate(sorted(group['guideID'].unique()), start=0)}
    
    # Apply the mapping
    group['enhancerID'] = group['enhancerID'].map(enhancer_mapping)
    group['guideID'] = group['guideID'].map(guide_mapping)
    
    print(group['enhancerID'].unique())
    print(group['guideID'].unique())
    
    # Prepare data for Stan model
    stan_data = {
        'N': len(group),
        'E': group['enhancerID'].nunique() - 1,  # Excluding control (-1)
        'G': group['guideID'].nunique() - 1,     # Excluding control (-1)
        'expression': group['expression'].astype(float).tolist(),
        'enhancer_idx': group['enhancerID'].tolist(),
        'guide_idx': group['guideID'].tolist()
    }

    # Save the data to a JSON file
    with open(output_file, 'w') as f:
        json.dump(stan_data, f)
