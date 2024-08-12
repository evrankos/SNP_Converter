from cli.arg_parser import parse_arguments
import pandas as pd
import numpy as np

def vcf_to_plink(input_file1, output_file1, output_file2):
    vcf_file = input_file1

    # Initialize an empty list to store all lines of data
    vcf_data = []
    vcf_columns = []

    with open(vcf_file, 'r') as f:
        first_line = True
        for line in f:
            if line.startswith('##'):
                continue  # Skip lines starting with '##'
            fields = line.strip().split('\t')
            if first_line:
                first_line = False
                vcf_columns = [col.lstrip('#') for col in fields]  # Remove '#' from column headers
                continue  # Skip the first line as it's the header row
            vcf_data.append(fields)  # Append the fields as a list to vcf_data

    # Create a DataFrame from the vcf_data with proper column names
    vcf_df = pd.DataFrame(vcf_data, columns=vcf_columns)
    
    verify_col = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    
    if verify_col == vcf_columns[:9]:
        pass
    else:
        return "Invalid vcf file. Check if the file is formatted correctly."
    
    # Create the map_df DataFrame
    map_df = pd.DataFrame({
        'CHROM': vcf_df['CHROM'],
        'ID': vcf_df['ID'],
        'CM': -9,
        'POS': vcf_df['POS']
    })
    
    # Define the path where you want to save the .plk.map file
    output_2 = output_file2

    # Save the map_df to a .plk.map file
    map_df.to_csv(output_2, sep='\t', header=False, index=False)
    
    df = vcf_df
    
    format = vcf_df.iloc[0, 8]

    if 'GT' in format:
        if ':' in format:
            gt_index = 0
        else:
            pass
    else:
        return "GenoType data not found."
    
    columns_to_modify = vcf_df.columns[2:]

    if gt_index == 0:
        # Apply the transformation
        for col in columns_to_modify:
            vcf_df[col] = vcf_df[col].apply(lambda x: x.split(':')[0] if len(x.split(':')) > 0 else x)
    
    def replace_value(row, sample_col):
        ref = row['REF']
        alt = row['ALT']
        value = row[sample_col]
        
        if value in ['./.', '.|.']:
            return '00'
        
        if value in ['0/0', '0|0']:
            return ref + ref
        elif value in ['1/1', '1|1']:
            return alt + alt
        elif value in ['1/0', '1|0']:
            return alt + ref
        elif value in ['0/1', '0|1']:
            return ref + alt
            
        else:
            return '00'
    
    # Apply replacement function to each sample column
    sample_columns = df.columns[9:]
    for col in sample_columns:
        df[col] = df.apply(lambda row: replace_value(row, col), axis = 1)
        
    df.drop(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)

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
    
    print("Your Variant Call Format file (.vcf) has been converted to PLINK files (.plk.ped and .plk.map).")
