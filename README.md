# SNP Converter

This Python package provides a command-line tool for converting genotype data files between PED, VCF, and HMP formats.

## Directory Structure

```plaintext
SNP_Converter/
├── cli/
│   ├── __init__.py
│   ├── arg_parser.py
├── conversions/
│   ├── __init__.py
│   ├── plink_to_hapmap.py
│   ├── plink_to_vcf.py
│   ├── hapmap_to_plink.py
│   ├── hapmap_to_vcf.py
│   ├── vcf_to_plink.py
│   ├── vcf_to_hapmap.py
├── utils/
│   ├── __init__.py
│   ├── file_parsers.py
├── main.py
├── README.md
├── requirements.txt
