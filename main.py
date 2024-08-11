from cli.arg_parser import parse_arguments
from conversions.plink_to_hapmap import plink_to_hmp
from conversions.plink_to_vcf import plink_to_vcf
from conversions.hapmap_to_plink import hmp_to_plink
from conversions.hapmap_to_vcf import hmp_to_vcf
from conversions.vcf_to_plink import vcf_to_plink
from conversions.vcf_to_hapmap import vcf_to_hmp

def main():
    args = parse_arguments()

    if args.operation == "plink_to_vcf":
        plink_to_vcf(args.input_file1, args.input_file2, args.output_file1)
    elif args.operation == "plink_to_hmp":
        plink_to_hmp(args.input_file1, args.input_file2, args.output_file1)
    elif args.operation == "hmp_to_plink":
        hmp_to_plink(args.input_file1, args.output_file1, args.output_file2)
    elif args.operation == "hmp_to_vcf":
        hmp_to_vcf(args.input_file1, args.output_file1)
    elif args.operation == "vcf_to_plink":
        vcf_to_plink(args.input_file1, args.output_file1, args.output_file2)
    elif args.operation == "vcf_to_hmp":
        vcf_to_hmp(args.input_file1, args.output_file1)
    
if __name__ == "__main__":
    main()
