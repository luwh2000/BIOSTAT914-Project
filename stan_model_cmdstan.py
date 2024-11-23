import pandas as pd
import os

def parse_stan_output(file_path):
    """
    Parses the Stan output CSV file and returns a Pandas DataFrame.
    
    Parameters:
        file_path (str): Path to the Stan output CSV file.
    
    Returns:
        pd.DataFrame: DataFrame containing the sampled values.
    """
    # Open the file and read lines
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Extract the header and data start line
    header_line = None
    data_start_index = None
    
    for i, line in enumerate(lines):
        if line.startswith("#"):
            continue  # Skip comments
        if header_line is None:
            header_line = line.strip()  # First non-comment line is the header
        else:
            data_start_index = i
            break

    # Remove lines that start with '#'
    lines = [line for line in lines if not line.startswith("#")]
    
    # Ensure header and data were found
    if header_line is None or data_start_index is None:
        raise ValueError("Invalid Stan output file format.")
    
    # Read the header and the data
    headers = header_line.split(",")
    data_lines = lines[data_start_index:]
    # print(data_lines)
    
    # Convert the data into a list of lists
    data = [list(map(float, line.strip().split(","))) for line in data_lines if line.strip()]
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=headers)
    
    return df


# Example usage
if __name__ == "__main__":
    # Replace 'output/results.csv' with the actual path to your Stan output file
    output_folder = "output"

    for i in range(1, 2):
        out_file = f"output_{i}.csv"
        stan_output_file = os.path.join(output_folder, output_file)
        
        stan_df = parse_stan_output(stan_output_file)
        print(f"Parsed Stan Output from {output_file}:")
        print(stan_df.head())
    # average each column
    print(stan_df.mean())
