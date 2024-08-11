from cli.arg_parser import parse_arguments
import pandas as pd
import numpy as np

def plink_to_hapmap(input_file1, input_file2, output_file1):
    # Define the path to .plk.map file
    map_file = input_file2

    # Read the .plk.map file into a dataframe
    map_df = pd.read_csv(map_file, sep='\t', header=None, names=['CHROM', 'ID', 'CM', 'POS'])
    
    # Define the path to .plk.ped file
    ped_file = input_file1
    
    # Read the .plk.ped file into a dataframe
    ped_df = pd.read_csv(ped_file, sep='\t', header=None)
    
    df = ped_df.copy()
    
    hmp_df = df.T
    
    # Drop specific rows by index
    rows_to_drop = [0, 2, 3, 4, 5]
    hmp_df = hmp_df.drop(rows_to_drop)

    hmp_df.columns = hmp_df.iloc[0]
    hmp_df = hmp_df[1:].reset_index(drop=True)

    # Initialize lists to store REF and ALT values
    ref_alleles = []
    alt_alleles = []

    # Process every pair of rows
    for i in range(0,  hmp_df.shape[0], 2):
        # Select the pair of rows
        row_pair =  hmp_df.iloc[i:i+2]
        
        # Flatten the values of both rows into a single list
        values = row_pair.values.flatten()
        
        # Filter out '0' values
        filtered_values = [val for val in values if val != '0' and pd.notna(val)]
        
        # Count frequency of each allele
        allele_counts = pd.Series(filtered_values).value_counts()
        
        # Most common allele is REF
        if len(allele_counts) > 0:
            ref_allele = allele_counts.idxmax()
            ref_alleles.append(ref_allele)
            
            # Second most common allele is ALT
            if len(allele_counts) > 1:
                alt_allele = allele_counts.drop(ref_allele).idxmax()
                alt_alleles.append(alt_allele)
            else:
                alt_alleles.append(np.nan)  # In case there's no ALT
        else:
            ref_alleles.append(np.nan)
            alt_alleles.append(np.nan)

    # Create a DataFrame with REF and ALT as the first two columns
    ref_alt_df = pd.DataFrame({'REF': ref_alleles, 'ALT': alt_alleles})

    # Concatenate every 2 rows into a single row
    df_combined =  hmp_df.groupby(np.arange(len(hmp_df)) // 2).apply(lambda x: x[::-1].astype(str).sum()).reset_index(drop=True)

    # Combine the REF/ALT DataFrame with the combined rows
    df_final = pd.concat([ref_alt_df, df_combined], axis=1)

    hmp_df = df_final

    # List of new columns to be added
    new_columns = ['strand', 'assembly#', 'center', 'protLSID', 'assayLSID', 'panel', 'QCcode']

    # Create a DataFrame with these columns and all values as NaN
    new_columns_df = pd.DataFrame(columns=new_columns)
    new_columns_df = new_columns_df.reindex(hmp_df.index)  # Make sure it has the same index as hmp_df

    # Concatenate the new columns DataFrame with hmp_df
    hmp_df = pd.concat([new_columns_df, hmp_df], axis=1)

    # Ensure the lengths of map_df and hmp_df match
    if len(map_df) != len(hmp_df):
        raise ValueError("map_df and hmp_df must have the same number of rows")

    # Add the 'pos' column to hmp_df
    hmp_df.insert(0, 'pos', map_df['POS'])

    # Add the 'chrom' column to hmp_df
    hmp_df.insert(0, 'chrom', map_df['CHROM'])

    hmp_df['alleles'] = hmp_df['ALT'] + '/' + hmp_df['REF']

    # Remove the 'REF' and 'ALT' columns
    hmp_df = hmp_df.drop(columns=['REF', 'ALT'])

    # Reorder columns to make 'alleles' the first column
    # Create a list of columns in the new order
    columns = ['alleles'] + [col for col in hmp_df.columns if col != 'alleles']

    # Reorder the DataFrame columns
    hmp_df = hmp_df[columns]

    # Add the 'ID' column to hmp_df
    hmp_df.insert(0, 'rs#', map_df['ID'])

    # Identify columns starting from the 11th column (0-based index 10)
    start_col_index = 11
    columns_to_modify = hmp_df.columns[start_col_index:]

    # Replace '00' with 'NN' in these columns
    hmp_df[columns_to_modify] = hmp_df[columns_to_modify].replace('00', 'NN')

    hmp_df = hmp_df.fillna('NA')

    # Specify the file path where you want to save the file
    file_path = output_file1

    # Write the DataFrame to a .hmp.txt file
    hmp_df.to_csv(file_path, sep='\t', index=False)

    print("Your PLINK files (.plk.ped and .plk.map) has been converted to a HapMap file (.hmp.txt).")
