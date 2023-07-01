#!/usr/bin/env python3

# *-----------------------------------------------------------------------------------------------
# | PROGRAM NAME: extract_from_genbank.py
# | DATE: 26/08/22 
# | CREATED BY: Lila Maciel Rodriguez Perez
# | PROJECT FILE: get_from_genbank
# *-----------------------------------------------------------------------------------------------
# | PURPOSE: Extract CDS dnam CDS protein, rRNA or tRNA from genbank files annotated by Prokka
# *------------------------------------------------------------------------------------------
# | USAGE:  extract_from_genbank.py -ft 'CDS' -name 'pbp2' -i data/*.gbk -o outputs/pbp2.fasta 
# |         extract_from_genbank.py -ft 'protein' -name 'pbp2' -i data/*.gbk -o outputs/pbp2.fasta 
# *-----------------------------------------------------------------------------------------------
# | WARNING: Output directory must be inside the directory in which we are running the script
# *-----------------------------------------------------------------------------------------------

# *-------------------------------------  Libraries -----------------------------------------------
import os # SyntaxError: import * only allowed at module level
import sys
import argparse

#from pyfiglet import Figlet

from Bio import SeqIO
from Bio.Seq import * # funcions like reverse_complement, translate, etc
from Bio.SeqRecord import *  ####====> to get sequence in fasta format


# *--------------------------------------- Parsing arguments --------------------------------------
# creating parser. Adding description
parser = argparse.ArgumentParser(description='Program to extract any feature from 1 or more genbank files.')
# Adding arguments
parser.add_argument('-ft', '--feature_type', type = str, help = 'Can be gene (CDS dna seq), protein, rRNA or tRNA')
parser.add_argument('-n', '--name', type = str, help = 'Name of the gene or the protein product or symbol.')
# input_gb is a list
parser.add_argument('-i', '--input_gb', 
                    nargs = '+', 
                    type = str,
                    help='Genbank files from which to extract the specified sequences')
parser.add_argument('-o', '--output', 
                    type=str, 
                    help='File in which user wants to store the results')
# Parsing arguments
args = parser.parse_args()

# *--------------------------------------- Defining functions -------------------------------------

def get_file_basename(input_file):
    """
    Returns the name of the file
    without the path of directories to the file
    """
    # variables
    file_path = input_file
    file_basename = os.path.basename(input_file)
    file_path_dirs = os.path.dirname(file_path)
    file_dirname = os.path.basename(file_path_dirs)
    return file_basename

def get_directory_name(input_file):
    """
    Returns the name of the directory where the input_file is.
    """
    # variables
    file_path = input_file
    file_basename = os.path.basename(input_file)
    file_path_dirs = os.path.dirname(file_path)
    file_dirname = os.path.basename(file_path_dirs)
    return file_dirname

def dna_or_protein(answer, sequence):
    """
    Test if user wants dna or protein sequence
    
    Parameters: 
    * answer to dna or protein?: Boolean True or False
    * sequence: dna sequence
    
    Returns dna or protein sequence
    """ 
    
    if answer == False:
        # User wants dna sequence
        #sequence = dna_sequence
        return sequence
                    
    else:
        # User wants protein sequence
        # Translate the dna sequence
        protein_sequence = translate(
            sequence=sequence,
            #to_stop=True,
            #cds=True,
            to_stop=False,
            table=11
            )
        return protein_sequence
    
def extract_feature_from_genbank(feature_type,name,input_gb): #,output_fasta):
    """
    Extract dna or protein sequence of CDS, rRNA, tRNA
    
    Parameters: 
    * feature_type: CDS (dna), protein, rRNA, tRNA.
    * name: of gene, protein product, rRNA or tRNA.
    * input_gb: input genbank file.
    * Result: sequence list of a single genbank file
    """
    # Variables
    ## List to store SeqRecord objects (sequence) of a single genbank file
    sequence_list = []
    ## flag to test if user wants gene or protein
    protein_wanted_flag = False
    ## test if user want dna or protein sequence:
    if feature_type == 'protein':
        print("User wants protein!")
        protein_wanted_flag = True
        feature_type = 'CDS'
        # IN THIS CASE, we are always going to look in CDS
    
    # Start of search processs
    record = SeqIO.parse(input_gb, "genbank")

    for element in record:
        for feature in element.features:
            
            ftype = feature.type
            # IN THIS CASE, we are always going to look in CDS
            # CDS has a quealifiers dictionary with 'product' as key, generally
            
            if 'product' not in feature.qualifiers.keys():
                continue # to the next feature iteration
            else:
                product = feature.qualifiers['product'][0] ### LIST, RIGHT??
            
            if not ( ftype == feature_type and ( product == name or product == f"{name} (partial)" ) ):
                continue # to the next feature iteration
            else:
                # if the type and name of product are both what I wanted, start processing:
                contig_id = element.id
                location = str(feature.location)
                start = feature.location.start
                end = feature.location.end
                strand = feature.location.strand
                # specially for mishell directories ===> change later
                sample = get_directory_name(input_gb) 
                
                if strand == +1:
                    dna_seq = element.seq[start:end] # CDS
                    
                    # Get the dna or protein sequence with dna_or_protein()
                    new_seq = dna_or_protein(answer=protein_wanted_flag, sequence=dna_seq)
                    len_seq = len(new_seq)
                    # Generate SeqRecord object
                    seq_record = SeqRecord(
                        seq=new_seq,
                        id=f"Sample: {sample} | Product: {name} | Length: {len_seq} |",
                        description=f"Location: {location} | Contig: {contig_id}"
                        )
                    # Append the NEW SeqRecord object to list
                    sequence_list.append(seq_record)

                else: # -1
                    dna_seq = reverse_complement(element.seq[start:end]) # CDS
                    
                    # Get the dna or protein sequence with dna_or_protein()
                    new_seq = dna_or_protein(answer=protein_wanted_flag, sequence=dna_seq)
                    len_seq = len(new_seq)
                    # Generate SeqRecord object to store the sequence
                    seq_record = SeqRecord(
                        seq=new_seq,
                        id=f"Sample: {sample} | Product: {name} | Length: {len_seq} |",
                        description=f"Location: {location} | Contig: {contig_id}"
                        )
                    # Append the NEW SeqRecord object to list
                    sequence_list.append(seq_record)
    
    return sequence_list # of a single genbank file

# *------------------------------------------ Main script -----------------------------------------

# Arguments
feature_type = args.feature_type
name = args.name
## args.input_gb is a list of strings
genbank_file_list = args.input_gb # each element (input_gb) is a string
results_file = args.output


# variables
required_feature_list = []
counter = 0

# Process
for gb_file in genbank_file_list:
    counter += 1
    sample = get_directory_name(gb_file)
    
    print("\n")
    print(f"Processing sample #{counter}: {sample}")
    result_list = extract_feature_from_genbank(feature_type,name,gb_file)
    
    # whas the feature found in the current file?
    if len(result_list) < 1:
        print(f"{name} was not found in {sample}.")
    else:
        print(f"{name} was found in {sample}.")
        required_feature_list.extend(result_list)
    
    print("\n")
    # to the next genbank file

print(f"Total number of genbank files analized: {counter}")
print(f"Number of {name} sequences found in all the gb files: {len(required_feature_list)}")
print("\n")

SeqIO.write(required_feature_list, results_file, "fasta")

print(f"Finished!!! Happy analysis of your {name} sequences!!")
print("\n")
