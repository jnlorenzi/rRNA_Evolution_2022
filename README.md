# rRNA_Evolution_2022
Script used for the rRNA evolution paper (2022)

## requirements

### _Python_ mods
os

sys

argparse

json



## script description

for each python script, the command : python [script_name].py --help returns detailed execution information

_make_core.py_ : from pairwise blastp, returns homolog and core genome relationship.

## pipeline

1.
manual recuperation of the genome fna file and annotation gbff file on the NCBI website, see report.txt file respectively in pseudomonas_2022/data/genome and pseudomonas_2022/data/annotation for more information.

2.
python annotation_and_genome_files_match.py /perso/lorenzi/Documents/pseudomonas_2022/data/genome/ /perso/lorenzi/Documents/pseudomonas_2022/data/annotation/ /perso/lorenzi/Documents/pseudomonas_2022/data/match_species.tab -g fna -a gbff

3.
python genbank_to_json.py /perso/lorenzi/Documents/pseudomonas_2022/data/annotation/ /perso/lorenzi/Documents/pseudomonas_2022/data/json/ -a gbff -i /perso/lorenzi/Documents/pseudomonas_2022/data/match_species.tab

4.
python species_filter.py /perso/lorenzi/Documents/pseudomonas_2022/data/json/ /perso/lorenzi/Documents/pseudomonas_2022/data/result/ -j txt






## archive 
python retrieve_genome_NCBI.py ~/edirect/ /perso/lorenzi/Documents/pseudomonas_2022/data/assembly_result.txt /perso/lorenzi/Documents/pseudomonas_2022/data/genome/ fasta
python retrieve_genome_NCBI.py ~/edirect/ /perso/lorenzi/Documents/pseudomonas_2022/data/assembly_result.txt /perso/lorenzi/Documents/pseudomonas_2022/data/genome/ gb


