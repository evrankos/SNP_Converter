import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert Genotype data files between PED, VCF, and HMP formats.")
    
    parser.add_argument("operation", choices=[
        "plink_to_vcf", 
        "plink_to_hapmap", 
        "vcf_to_plink", 
        "hapmap_to_plink",
        "hapmap_to_vcf",
        "vcf_to_hapmap"
    ], help="The conversion operation to perform ('plink_to_vcf', 'plink_to_hapmap', 'vcf_to_plink', 'hapmap_to_plink', 'hapmap_to_vcf', 'vcf_to_hapmap').")
    
    parser.add_argument("input_file1", help="Path to the first input SNP file")
    parser.add_argument("input_file2", nargs='?', help="Path to the second input SNP file (required for specific operations)")
    parser.add_argument("output_file1", help="Path to the first output SNP file")
    parser.add_argument("output_file2", nargs='?', help="Path to the second output SNP file (produced for specific operations)")
    
    args = parser.parse_args()

    # Validate arguments based on the operation
    if args.operation in ["plink_to_vcf", "plink_to_hapmap"]:
        if not args.input_file2:
            parser.error(f"Operation '{args.operation}' requires two input files (PED and MAP).")
        if args.output_file2:
            parser.error(f"Operation '{args.operation}' produces only one output file, do not specify the second output file.")
    
    if args.operation in ["vcf_to_plink", "hapmap_to_plink"]:
        if not args.output_file2:
            parser.error(f"Operation '{args.operation}' requires two output files (PED and MAP).")
        if args.input_file2:
            parser.error(f"Operation '{args.operation}' requires only one input file, do not specify the second input file.")
            
    if args.operation in ["plink_to_vcf", "plink_to_hapmap", "vcf_to_plink", "vcf_to_hapmap", "hapmap_to_plink", "hapmap_to_vcf"]:
        if not args.input_file1:
            parser.error(f"Input file (path) is missing.")
        if not args.output_file1:
            parser.error(f"Output file (path) is missing.")

    return args
