from cli.arg_parser import parse_arguments
import pandas as pd
import numpy as np

def plink_to_vcf(input_file1, input_file2, output_file1):
    # Define the path to .plk.map file
    map_file = input_file2

    # Read the .plk.map file into a dataframe
    map_df = pd.read_csv(map_file, sep='\t', header=None, names=['CHROM', 'ID', 'CM', 'POS'])
    
    # Define the path to .plk.ped file
    ped_file = input_file1
    
    # Read the .plk.ped file into a dataframe
    ped_df = pd.read_csv(ped_file, sep='\t', header=None)
    
    df = ped_df.copy()
    
    vcf_df = df.T
    
    # Drop specific rows by index
    rows_to_drop = [0, 2, 3, 4, 5]
    vcf_df = vcf_df.drop(rows_to_drop)
    
    vcf_df.columns = vcf_df.iloc[0]
    vcf_df = vcf_df[1:].reset_index(drop=True)

    # Initialize lists to store REF and ALT values
    ref_alleles = []
    alt_alleles = []

    # Process every pair of rows
    for i in range(0,  vcf_df.shape[0], 2):
        # Select the pair of rows
        row_pair =  vcf_df.iloc[i:i+2]
        
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
                alt_alleles.append('.')  # In case there's no ALT
        else:
            ref_alleles.append(np.nan)
            alt_alleles.append(np.nan)

    # Create a DataFrame with REF and ALT as the first two columns
    ref_alt_df = pd.DataFrame({'REF': ref_alleles, 'ALT': alt_alleles})

    # Concatenate every 2 rows into a single row
    df_combined =  vcf_df.groupby(np.arange(len(vcf_df)) // 2).apply(lambda x: x[::-1].astype(str).sum()).reset_index(drop=True)

    vcf_df = df_combined

    # List of new columns to be added
    new_columns = ["QUAL", "FILTER", "INFO", "FORMAT"]

    # Create a DataFrame with these columns and all values as NaN
    new_columns_df = pd.DataFrame(columns=new_columns)
    new_columns_df = new_columns_df.reindex(vcf_df.index)  # Make sure it has the same index as vcf_df

    # Concatenate the new columns DataFrame and REF/ALT DataFrame with vcf_df
    vcf_df = pd.concat([ref_alt_df, new_columns_df, vcf_df], axis=1)

    # Add the 'ID' column to vcf_df
    vcf_df.insert(0, 'ID', map_df['ID'])

    # Add the 'pos' column to vcf_df
    vcf_df.insert(0, 'POS', map_df['POS'])

    # Add the 'chrom' column to vcf_df
    vcf_df.insert(0, 'CHROM', map_df['CHROM'])

    # Identify columns starting from the 11th column (0-based index 10)
    start_col_index = 9
    columns_to_modify = vcf_df.columns[start_col_index:]

    # Replace '00' with 'NN' in these columns
    vcf_df[columns_to_modify] = vcf_df[columns_to_modify].replace('00', 'NN')

    vcf_df['FORMAT'] = 'GT'
    vcf_df['FILTER'] = 'PASS'

    def replace_value(row, sample_col):
        ref = row['REF']
        alt = row['ALT']
        value = row[sample_col]
        
        if value == 'NN':
            return './.'
            
        allele1 = value[0]
        allele2 = value[1] 
        
        if allele1 == allele2:
            if allele1 == ref:
                return '0/0'
            elif allele1 == alt:
                return '1/1'
            else:
                return './.'
        else:
            if allele1 == ref and allele2 == alt:
                return '1/0'
            elif allele1 == alt and allele2 == ref:
                return '0/1'
            else:
                return './.'
            
    # Apply replacement function to each sample column
    sample_columns = vcf_df.columns[9:]
    for col in sample_columns:
        vcf_df[col] = vcf_df.apply(lambda row: replace_value(row, col), axis = 1)
        
    # The VCF standard does not permit the REF allele to be missing. As a consequence, PLINK converts missing REF alleles (which can appear in e.g. data imported from PLINK 1 .ped files) to 'N' when exporting VCF files.
    vcf_df['REF'] = vcf_df['REF'].fillna('N') 

    vcf_df = vcf_df.fillna('.')

    # Define the output file path
    output_file = output_file1

    # Define the header for the VCF file
    header = """##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"""

    # Prepare the header by adding sample names from the DataFrame
    header += '\t' + '\t'.join(vcf_df.columns[9:]) + '\n'

    # Function to convert a row of DataFrame to VCF format
    def row_to_vcf(row):
        chrom = str(row['CHROM'])
        pos = str(row['POS'])
        _id = row['ID']
        ref = row['REF']
        alt = row['ALT']
        qual = row['QUAL']
        filt = row['FILTER']
        info = row['INFO']
        format_field = 'GT'
        genotypes = '\t'.join(row.iloc[9:].tolist())
        
        return f"{chrom}\t{pos}\t{_id}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{format_field}\t{genotypes}"

    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Write the header
        f.write(header)
        
        # Write each row in VCF format
        for index, row in vcf_df.iterrows():
            vcf_line = row_to_vcf(row)
            f.write(vcf_line + '\n')
            
    print("Your PLINK files (.plk.ped and .plk.map) has been converted to a Variant Call Format file (.vcf).")