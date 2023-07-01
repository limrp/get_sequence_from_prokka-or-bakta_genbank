# get_sequence_from_prokka-or-bakta_genbank

This script allows you to extract gene and protein sequences from one or more GenBank files.

## Dependencies

This script requires Python 3 and the following Python libraries:

- Biopython 1.79
- os 
- sys
- argparse

## Command Line Arguments

The script accepts the following command line arguments:

- `-h`, `--help`: Show the help message and exit.
- `-ft FEATURE_TYPE`, `--feature_type FEATURE_TYPE`: The type of feature to extract. Can be 'gene' (CDS DNA sequence), 'protein', 'rRNA', or 'tRNA'.
- `-n NAME`, `--name NAME`: The name of the gene or the protein product or symbol.
- `-i INPUT_GB [INPUT_GB ...]`, `--input_gb INPUT_GB [INPUT_GB ...]`: The GenBank files from which to extract the specified sequences.
- `-o OUTPUT`, `--output OUTPUT`: The file in which to store the results.

## Usage

You can run the script using the following command:

```bash
python3 get_seq_from_genbank.py -ft 'CDS' -n 'cell division protein (ftsH)' -i test_data/*/*.gbff -o ftsH_test_results_CDS.fasta
