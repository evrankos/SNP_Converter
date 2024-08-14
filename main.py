from cli.arg_parser import parse_arguments
from conversions.plink_to_hapmap import plink_to_hapmap
from conversions.plink_to_vcf import plink_to_vcf
from conversions.hapmap_to_plink import hapmap_to_plink
from conversions.hapmap_to_vcf import hapmap_to_vcf
from conversions.vcf_to_plink import vcf_to_plink
from conversions.vcf_to_hapmap import vcf_to_hapmap

def main():
    args = parse_arguments()

    if args.operation == "plink_to_vcf":
        plink_to_vcf(args.input_file1, args.input_file2, args.output_file1)
    elif args.operation == "plink_to_hapmap":
        plink_to_hapmap(args.input_file1, args.input_file2, args.output_file1)
    elif args.operation == "hapmap_to_plink":
        hapmap_to_plink(args.input_file1, args.output_file1, args.output_file2)
    elif args.operation == "hapmap_to_vcf":
        hapmap_to_vcf(args.input_file1, args.output_file1)
    elif args.operation == "vcf_to_plink":
        vcf_to_plink(args.input_file1, args.output_file1, args.output_file2)
    elif args.operation == "vcf_to_hapmap":
        vcf_to_hapmap(args.input_file1, args.output_file1)
    
if __name__ == "__main__":
    main()
