from cli.arg_parser import parse_arguments
import pandas as pd

def hapmap_to_vcf(input_file1, output_file1):

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
    
    verify_col = ["rs", "alleles", "chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panel", "QCcode"]
    
    if verify_col == hmp_columns[:11]:
        pass
    else:
        return "Invalid hmp.txt file. Check if the file is formatted correctly."
    
    vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    for i in range(11, len(hmp_columns)):
        if hmp_columns[i] not in vcf_columns:
            vcf_columns.append(hmp_columns[i])     

    vcf_df = pd.DataFrame(index=range(len(hmp_df)), columns=vcf_columns)

    # Extract the Sample ID columns from vcf_df
    hmp_sample_ids = hmp_df.iloc[:, 11:]

    # Update the Sample ID columns of hmp_df with the values from vcf_df
    vcf_df.iloc[:, 9:] = hmp_sample_ids.values

    vcf_df['ID'] = hmp_df['rs']
    vcf_df['ALT'] = hmp_df['alleles'].str.split('/').str[0]
    vcf_df['REF'] = hmp_df['alleles'].str.split('/').str[1]
    vcf_df['CHROM'] = hmp_df['chrom']
    vcf_df['POS'] = hmp_df['pos']

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
            
    print("Your HapMap file (.hmp.txt) has been converted to a Variant Call Format file (.vcf).")
