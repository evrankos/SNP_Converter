from cli.arg_parser import parse_arguments
import pandas as pd

def vcf_to_hmp(input_file1, output_file1):
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

    hmp_columns = ["rs#", "alleles", "chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panel", "QCcode"]

    for i in range(9, len(vcf_columns)):
        if vcf_columns[i] not in hmp_columns:
            hmp_columns.append(vcf_columns[i])     

    hmp_df = pd.DataFrame(index=range(len(vcf_df)), columns=hmp_columns)

    # Extract the Sample ID columns from vcf_df
    vcf_sample_ids = vcf_df.iloc[:, 9:]

    # Update the Sample ID columns of hmp_df with the values from vcf_df
    hmp_df.iloc[:, 11:] = vcf_sample_ids.values

    hmp_df['rs#'] = vcf_df['ID']
    hmp_df['alleles'] = vcf_df['ALT'] + '/' + vcf_df['REF']
    hmp_df['chrom'] = vcf_df['CHROM']
    hmp_df['pos'] = vcf_df['POS']

    def replace_values(row, alleles):
        # Precompute allele mappings
        allele1, allele2 = alleles.split('/')
        allele1_parts = allele1.split(',')
        
        if len(allele1_parts) > 1:
            allele1 = allele1_parts[0]
            allele12 = allele1_parts[1]
            allele13 = allele1_parts[2]
            
            allele_map = {
            '0/0': f"{allele2}{allele2}",
            '1/1': f"{allele1}{allele1}",
            './.': 'NN',
            '1/0': f"{allele2}{allele1}",
            '0/1': f"{allele1}{allele2}",
            '2/2': f"{allele12}{allele12}",
            '3/3': f"{allele13}{allele13}",
            '3/2': f"{allele12}{allele13}",
            '2/3': f"{allele13}{allele12}"
        }
        else:
            allele_map = {
            '0/0': f"{allele2}{allele2}",
            '1/1': f"{allele1}{allele1}",
            './.': 'NN',
            '1/0': f"{allele2}{allele1}",
            '0/1': f"{allele1}{allele2}"
        }
        
        # Use row.map() with a default value 'NN' for any values not in allele_map
        return row.map(lambda x: allele_map.get(x, 'NN'))

    for index, row in hmp_df.iterrows():
        alleles = row['alleles']
        hmp_df.loc[index, hmp_df.columns[11:]] = replace_values(row[hmp_df.columns[11:]], alleles)
        
    hmp_df = hmp_df.fillna('NA')

    # Specify the file path where you want to save the file
    file_path = output_file1

    # Write the DataFrame to a .hmp.txt file
    hmp_df.to_csv(file_path, sep='\t', index=False)
    
    print("Your Variant Call Format file (.vcf) has been converted to a HapMap file (.hmp.txt).")
