import pandas as pd

# Step 1: Load the data
file_path = "data1.txt"  # Replace with your actual file path
data = pd.read_csv(file_path, sep="\t")

data = data[(data['geneID'] == 1) & (data['enhancerID'] == 1) & (data['guideID'] == 1)]
print(data)
