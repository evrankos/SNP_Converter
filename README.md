# SNP Converter

## Overview

The SNP Converter package provides a command-line tool to convert genotype data files between three formats: PED+MAP, VCF, and HMP. This package is designed to assist geneticists, bioinformaticians, and researchers in performing these conversions efficiently and accurately.

## Directory Structure
The package is organized as follows:
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

## Components

- `cli/arg_parser.py`: Handles argument parsing for the command-line tool.
- `conversions/`: Contains scripts for specific file format conversions.
- `utils/file_parsers.py`: Placeholder for any file parsing utilities or helper functions.
- `main.py`: Contains the main function that calls the appropriate conversion function based on user-specified operations.

## Installation

#### Prerequisites

- **Python Version:** Ensure that Python 3.6 or later is installed on your system. You can download Python from the official website [here](https://www.python.org/downloads/).

### Step 1: Clone the Repository

First, clone the repository to your local machine:
```bash session
git clone https://github.com/evrankos/SNP_Converter.git
```

Navigate to the project directory:
```bash session
cd path/to/SNP_Converter
```

### Step 2: Create a Virtual Environment (Optional but Recommended)

It's recommended to create a virtual environment to avoid conflicts with other Python packages:
```bash session
python -m venv env
```

Activate the virtual environment:
- **Windows:**
  ```bash session
  .\env\Scripts\activate
  ```
- **macOS/Linux:**
  ```bash session
  source env/bin/activate
  ```

### Step 3: Install the Required Packages

Use the `requirements.txt` file to install all necessary dependencies:
```bash session
pip install -r requirements.txt
```

This command will install the following packages:
- `argparse==1.1`
- `numpy==2.0.1`
- `pandas==2.2.2`

### Step 4: Run the Command-Line Tool

After installing the dependencies, you can run the command-line tool provided in the package. Navigate to the package directory and execute the script with appropriate arguments.

For example, to convert a VCF file to a HapMap file, you can run:
```bash session
python main.py vcf_to_hmp /path/to/input.vcf /path/to/output.hmp.txt
```

## Additional Notes

- Make sure that your input files are properly formatted and located in accessible paths.
- You can find more detailed usage instructions in the "Usage" section of this documentation.

## Usage

You can use the SNP Converter tool via the command line. The main script is `main.py`, and it can be invoked with specific options for each type of conversion.

### Command-Line Interface

Use the following command to run the tool:
```bash session
python main.py <operation> <input_file> <output_file> [additional_options]
```
- `<operation>`: Specify the operation to perform. Options include:
  - `plink_to_hapmap`
  - `plink_to_vcf`
  - `hapmap_to_plink`
  - `hapmap_to_vcf`
  - `vcf_to_plink`
  - `vcf_to_hapmap`
- `<input_file>`: The input genotype file.
- `<output_file>`: The output file for the converted data.
- `[additional_options]`: Any additional options required by specific operations.

### Example Commands

1. **Convert PLINK to HapMap:**
   ```bash session
   python main.py plink_to_hapmap input.plk.ped input.plk.map output.hmp.txt
   ```
2. **Convert PLINK to VCF:**
   ```bash session
   python main.py plink_to_vcf input.plk.ped input.plk.map output.vcf
   ```
3. **Convert HapMap to PLINK:**
   ```bash session
   python main.py hapmap_to_plink input.hmp.txt output.plk.ped output.plk.map
   ```
4. **Convert HapMap to VCF:**
   ```bash session
   python main.py hapmap_to_vcf input.hmp.txt output.vcf
   ```
5. **Convert VCF to PLINK:**
   ```bash session
   python main.py vcf_to_plink input.vcf output.plk.ped output.plk.map
   ```
6. **Convert VCF to HapMap:**
   ```bash session
   python main.py vcf_to_hapmap input.vcf output.hmp.txt
   ```

## Conversion Scripts

### `plink_to_hapmap.py`

This script converts PLINK files into HapMap format.
- **`plink_to_hapmap(input_file1, input_file2, output_file1)`:** Converts PLINK files to a HapMap file.
  - **Input:**
    - **`input_file1`:** Path to the PLINK PED file.
    - **`input_file2`:** Path to the PLINK MAP file.
    - **`output_file1`:** Path to the output HapMap file.
  - **Process:**
    - Reads and processes the PLINK PED and MAP data.
    - Validates columns.
    - Converts genotype data.
    - Exports the data to a HapMap format file.
  - **Output:** A `.hmp.txt` file containing the converted genotype data.

### `plink_to_vcf.py`

This script converts PLINK files into VCF format.
- **`plink_to_vcf(input_file1, input_file2, output_file1)`:** Converts PLINK files to a VCF file.
  - **Input:**
    - **`input_file1`:** Path to the PLINK PED file.
    - **`input_file2`:** Path to the PLINK MAP file.
    - **`output_file1`:** Path to the output VCF file.
  - **Process:**
    - Reads and processes the PLINK PED and MAP data.
    - Validates columns.
    - Converts genotype data.
    - Exports the data to a VCF format file.
  - **Output:** A `.vcf` file containing the converted genotype data.

### `hapmap_to_plink.py`

This script converts HapMap file into PLINK format.
- **`hapmap_to_plink(input_file1, output_file1, output_file2)`:** Converts a HapMap file to PLINK files.
  - **Input:**
    - **`input_file1`:** Path to the HapMap file.
    - **`output_file1`:** Path to the output PLINK PED file.
    - **`output_file2`:** Path to the output PLINK MAP file.
  - **Process:**
    - Reads and processes the HapMap data.
    - Validates columns.
    - Converts genotype data.
    - Saves the data into `.plk.ped` and `.plk.map` files.
  - **Output:** `.plk.ped` and `.plk.map` files containing the converted genotype data.

### `hapmap_to_vcf.py`

This script converts a HapMap file into VCF format.
- **`hapmap_to_vcf(input_file1, output_file1)`:** Converts a HapMap file to a VCF file.
  - **Input:**
    - **`input_file1`:** Path to the HapMap file.
    - **`output_file1`:** Path to the output VCF file.
  - **Process:**
    - Reads and processes the HapMap data.
    - Validates columns.
    - Converts genotype data.
    - Exports the data to a VCF format file.
  - **Output:** A `.vcf` file containing the converted genotype data.

### `vcf_to_plink.py`

This script converts a VCF file into PLINK format.
- **`vcf_to_plink(input_file1, output_file1, output_file2)`:** Converts a VCF file to PLINK files.
  - **Input:**
    - **`input_file1`:** Path to the VCF file.
    - **`output_file1`:** Path to the output PLINK PED file.
    - **`output_file2`:** Path to the output PLINK MAP file.
  - **Process:**
    - Reads and processes the VCF data.
    - Validates columns.
    - Converts genotype data.
    - Saves the data into `.plk.ped` and `.plk.map` files.
  - **Output:** `.plk.ped` and `.plk.map` files containing the converted genotype data.

### `vcf_to_hapmap.py`

This script converts a VCF file into HapMap format.
- **`vcf_to_hapmap(input_file1, output_file1)`:** Converts PLINK files to a HapMap file.
  - **Input:**
    - **`input_file1`:** Path to the VCF file.
    - **`output_file1`:** Path to the output HapMap file.
  - **Process:**
    - Reads and processes the VCF data.
    - Validates columns.
    - Converts genotype data.
    - Exports the data to a HapMap format file.
  - **Output:** A `.hmp.txt` file containing the converted genotype data.

### `file_parsers.py`

This module is a placeholder for any file parsing utilities or helper functions that may be used across the conversion scripts. Currently, it is not populated with specific code but can be expanded as needed for more complex parsing operations.

## Error Handling

The conversion scripts include basic error handling for:
- Invalid file formats.
- Missing genotype data.
- Incorrect or unrecognized columns.
These errors will prompt the user to check the input files for correctness and format them accordingly.

## Contributions

Contributions to improve the SNP Converter tool are welcome. Please follow the standard fork-and-pull model, making sure to test your changes thoroughly. Here's how you can help:

1. Fork the repository.
2. Create a new branch (e.g., feature-branch).
3. Make your changes.
4. Test your changes thoroughly.
5. Commit and push your changes.
6. Submit a pull request.




This documentation should help you understand how to use the SNP Converter package effectively. If you encounter any issues or have further questions, feel free to reach out!