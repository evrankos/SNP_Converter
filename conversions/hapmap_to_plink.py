from cli.arg_parser import parse_arguments
import pandas as pd
import numpy as np

def hapmap_to_vcf(input_file1, output_file1, output_file2):

    hmp_file = input_file1

    # Initialize an empty list to store all lines of data
    hmp_data = []
    hmp_columns = []

    with open(hmp_file, 'r') as f:
        first_line = True  # Flag to skip the first line
        for line in f:
            if line.startswith('##'):
                continue  # Skip lines starting with '##'
            fields = line.strip().split('\t')
            if first_line:
                first_line = False
                hmp_columns = [col.rstrip('#') for col in fields]  # Remove '#' from column headers 
                continue  # Skip the first line as it's the header row
            hmp_data.append(fields)  # Append the fields as a list to hmp_data

    # Create a DataFrame from the list of lists (hmp_data)
    hmp_df = pd.DataFrame(hmp_data, columns=hmp_columns)
    
    # Create the map_df DataFrame
    map_df = pd.DataFrame({
        'CHROM': hmp_df['chrom'],
        'ID': hmp_df['rs'],
        'CM': -9,
        'POS': hmp_df['pos']
    })

    # Define the path where you want to save the .plk.map file
    output_2 = output_file2

    # Save the map_df to a .plk.map file
    map_df.to_csv(output_2, sep='\t', header=False, index=False)

    # Identify columns starting from the 11th column (0-based index 10)
    start_col_index = 11
    columns_to_modify = hmp_df.columns[start_col_index:]

    # Replace 'NN' with '00' in these columns
    hmp_df[columns_to_modify] = hmp_df[columns_to_modify].replace('NN', '00')

    hmp_df.drop(columns=['rs', 'alleles', 'chrom', 'pos', 'strand', 'assembly', 'center', 'protLSID', 'assayLSID', 'panel', 'QCcode'], inplace=True)

    df = hmp_df

    # Preallocate the new DataFrame with twice the number of rows
    num_rows, num_cols = df.shape
    new_data = np.empty((num_rows * 2, num_cols), dtype=df.dtypes[0].type)

    # Step 1: Extract the first and second characters using vectorized operations
    first_chars = df.applymap(lambda x: x[0] if len(x) == 2 else x)
    second_chars = df.applymap(lambda x: x[1] if len(x) == 2 else '')

    # Step 2: Fill the new_data array
    new_data[::2] = second_chars.values   # Fill even-indexed rows with second characters
    new_data[1::2] = first_chars.values   # Fill odd-indexed rows with first characters

    # Step 3: Convert the NumPy array back to a DataFrame
    df = pd.DataFrame(new_data, columns=df.columns)

    # Create a DataFrame with -9 values
    first_row = pd.DataFrame([[-9] * len(df.columns)], columns=df.columns)

    # Step 1: Convert column names into a new DataFrame
    column_names_df = pd.DataFrame([df.columns.tolist()], columns=df.columns)

    # Create a DataFrame with -9 values
    new_rows = pd.DataFrame([[-9] * len(df.columns)] * 4, columns=df.columns)

    # Step 2: Concatenate the new row DataFrame with the original DataFrame
    df = pd.concat([first_row, column_names_df, new_rows, df], ignore_index=True)

    ped_df = df.T

    # Reset the index of df and discard the old index
    ped_df = ped_df.reset_index(drop=True)

    # Define the path where you want to save the .plk.map file
    output_1 = output_file1

    # Save the df to a .plk.ped file
    ped_df.to_csv(output_1, sep='\t', header=False, index=False)
    
    print("Your HapMap file (.hmp.txt) has been converted to PLINK files (.plk.ped and .plk.map).")