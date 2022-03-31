#!/usr/bin/env python3

import os
import re
import sys
import shutil
import argparse


def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description = 'This script computes all the pairwise ANIb for each genome in the given repository (sequence must be in fasta format). \n\t # Command line example: python make_ani.py C:/Users/Lorenzi/Documents/Streptomyces/fasta/ C:/Users/Lorenzi/Documents/Streptomyces/anib/')

    # Positional mandatory arguments
    parser.add_argument("fasta_directory", help="absolute or relative path to the fasta directory", type=str)
    parser.add_argument("output_anib", help="absolute or relative path to the outfile anib directory", type=str)

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.5 - 31.03.2022')

    # Parse arguments
    args = parser.parse_args()

    return args


def ANIb_calculator(path_to_genome, path_to_output):
        """ Computes the ANIb matrix for genomes in the repository. 
        ANI method follows the basic algorithm:
                - Align the genome of organism 1 against that of organism 2, and identify
                the matching regions
                - Calculate the percentage nucleotide identity of the matching regions, as
                an average for all matching regions
                Methods differ on: (1) what alignment algorithm is used, and the choice of
                parameters (this affects the aligned region boundaries); (2) what the input
                is for alignment (typically either fragments of fixed size, or the most
                complete assembly available); (3) whether a reciprocal comparison is
                necessary or desirable.
        ANIb: uses BLASTN to align 1000nt fragments of the input sequences
        
        This script takes as main input a directory containing a set of
        correctly-formatted FASTA multiple sequence files. All sequences for a
        single organism should be contained in only one sequence file. The names of
        these files are used for identification, so it would be advisable to name
        them sensibly.
        """
        
        os.system('nice -n 19 ./average_nucleotide_identity.py -i ' + path_to_genome + ' -o ' + path_to_output + ' -m ANIb -v -f')


''' Calculates the ANIb distance matrix for all the species in the collection '''


def main():
  
  args = parseArguments()
  fasta_directory = args.fasta_directory
  output_anib =  args.output_anib
  
  ANIb_calculator(fasta_directory, output_anib)

if __name__ == "__main__":
    main()
