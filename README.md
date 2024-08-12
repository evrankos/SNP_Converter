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
```
## Contributing

```plaintext
Contributions are welcome! Here's how you can help:

1. Fork the repository.
2. Create a new branch (e.g., feature-branch).
3. Make your changes.
4. Test your changes thoroughly.
5. Commit and push your changes.
6. Submit a pull request.
```

## Code Style

```plaintext
- Follow PEP 8 guidelines.
- Use descriptive variable and function names.
- Document your code.
```