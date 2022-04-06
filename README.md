# rRNA_Evolution_2022
Script used for the rRNA evolution paper (2022)

## requirements

BLAST+ executable in the $PATH, or available on the command line

### _Python_ mods
os

sys

argparse

json

Biopython


## script description

For each python script, the command : python [script_name].py --help returns detailed execution information

_make_core_genome.py_ : from pairwise blastp, returns homolog and core genome relationship.

_average_nucleotide_identity.py_ : contains all the functions to performed the ANIb computation

_make_ani.py_ : calculates the ANIb distance matrix for all the species in the given directory

_make_persistence.py_ : calculates the persistence for each cds of each species in the given homolog dictionary (generates an intermediate json file containing the orthlog relationships)
