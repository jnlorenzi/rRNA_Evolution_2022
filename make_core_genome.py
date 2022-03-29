#!/usr/bin/env python3

#### MODS

import os
import sys
import json
import argparse

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description = 'This script detect the homolog and core genome from the res of a blastp analysis. The blastp output must be in -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps' format. Genome name and gene id must be unique. \n\t # Command line example: python make_core_genome.py C:/Users/Lorenzi/Documents/Streptomyces/Blast_output/ Streptomyces C:/Users/Lorenzi/Documents/Streptomyces/')

    # Positional mandatory arguments
    parser.add_argument("blast_output", help="absolute or relative path to the blast output directory", type=str)
    parser.add_argument("run", help="run name", type=str)
    parser.add_argument("output", help="absolute or relative path to the outfile directory (2 .json files created)", type=str)
    parser.add_argument('-i', '--identity_threshold', help="homolog identity threshold", type=float, default = 40)
    parser.add_argument('-e', '--evalue_threshold', help="homolog evalue threshold", type=float, default = 1e-10)
    parser.add_argument('-a', '--alignment_threshold', help="homolog alignment threshold", type=float, default = 70)
    parser.add_argument('-v', '--variation_threshold', help="homolog sequence length variation threshold", type=float, default = 999)


    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0 - 06.02.2022')

    # Parse arguments
    args = parser.parse_args()

    return args


def homolog_finder(path_blast_res, identity_threshold = 40, alignment_threshold = 70, evalue_threshold = 1e-10, sequence_length_variation_threshold = 999):
        """ Find all the homolog genes for each pair of organisms of the collection. BLASTP format 7 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' is required.
        A pair of homologous genes is defined as best hit between two genes which fit thresholds.
        Default thresholds are:
                > % identity >= 40%
                > % alignement >= 70% (with the smallest sequence)
                > E-value <= 1e-10
                > % sequence length variation = 999% (no treshold)
        """
        homolog = {}
        nb_file = sum(1 for f in os.listdir(path_blast_res) if f.endswith('.bl'))
        i = 0
        print('# Homolog dertermination...')
        for blast_res_file in os.listdir(path_blast_res):
                i += 1
                print('\r\t' + str(i) + '/' + str(nb_file) + ' .bl files processed' , end = '')
                if blast_res_file.endswith('.bl'):
                                ref, tar = blast_res_file.split('-vs-')
                                tar = tar[:-3]
                                if ref not in homolog:
                                        homolog[ref] = {}
                                if tar not in homolog[ref]:
                                        homolog[ref][tar] = {}
                                with open(path_blast_res + blast_res_file, 'r') as handler:
                                        #handler = open(path_blast_res + blast_res_file, 'r')
                                        line = handler.readline()
                                        no_hit = False
                                        while line:
                                                if line.startswith('# Query:'):
                                                        query = line[9:].split()[0]
                                                        line = handler.readline()
                                                        no_hit = False
                                                else:
                                                        line = handler.readline()
                                                        continue
                                                while line and line.startswith('#'):    # First hit expecting, if not hit, the next query is analysed
                                                        if line.startswith('# Query:'):
                                                                no_hit = True
                                                                break
                                                        line = handler.readline()
                                                if ref != tar:
                                                        if no_hit or not line:
                                                                homolog[ref][tar][query] = 'NA'
                                                                continue
                                                        line_handler = line.split()
                                                        query, subject, identity, alignment_length, evalue, query_length, subject_length = line_handler[0], line_handler[1], float(line_handler[2]), float(line_handler[3]), float(line_handler[10]), float(line_handler[12]), float(line_handler[13])
                                                        if abs(query_length - subject_length) <= min(query_length, subject_length) * (sequence_length_variation_threshold / 100):
                                                                alignement = (alignment_length / min(query_length, subject_length)) * 100
                                                                if evalue <= evalue_threshold and identity >= identity_threshold and alignement >= alignment_threshold:
                                                                        homolog[ref][tar][query] = subject
                                                                else:
                                                                        homolog[ref][tar][query] = 'NA'
                                                        else:
                                                                homolog[ref][tar][query] = 'NA'
                                                else:
                                                        homolog[ref][tar][query] = query
                                                        continue
        # manually add the species vs species blast
        for species in homolog:
          homolog[species][species] = {}
          tar = list(homolog[species].keys())[0]
          i = 0
          while tar == species:
            i += 1
            tar = list(homolog[species].keys())[i]
          for query in homolog[species][tar]:
            homolog[species][species][query] = query
          print(species, len(homolog[species][species]))
        print('\nDone\n')
        return homolog


def core_finder(homolog, condition, reference = '', path = ''):
        """ Returns the core genome with the annotation of the reference from the homolog dictionary as a list of annotation. By default, the core is defined for all organisms.
        
        dictionary structure:
        
                [organism name] = core list
                
                > organism name: name of organism in collection like "Streptomyces_coelicolor_A3(2)"
                > core list = list of ID 
                
        Examples:
                >>> core['Streptomyces_coelicolor_A3(2)']
                ['Sco1_4884', 'Sco1_5243', 'Sco1_3183', 'Sco1_5706', 'Sco1_1644', 'Sco1_1302', 'Sco1_2167', 'Sco1_2592', 'Sco1_1618', 'Sco1_1602', 'Sco1_4231', 'Sco1_4098', 'Sco1_5628', 'Sco1_1377', 'Sco1_5338', ...]
        """
        
        print('\n# Data format conversion...')
        ortholog_by_family_redundant = {}
        for ref in homolog:
                for tar in homolog[ref]:
                        if not ref in ortholog_by_family_redundant:
                                ortholog_by_family_redundant[ref] = {}
                        if ref != tar:
                                for CDS_ref in homolog[ref][tar]:
                                        candidat = homolog[ref][tar][CDS_ref]
                                        if candidat in homolog[tar][ref] and homolog[tar][ref][candidat] == CDS_ref:
                                                if CDS_ref not in ortholog_by_family_redundant[ref]:
                                                        ortholog_by_family_redundant[ref][CDS_ref] = []
                                                        ortholog_by_family_redundant[ref][CDS_ref].append(CDS_ref + '--:--' + ref)
                                                ortholog_by_family_redundant[ref][CDS_ref].append(candidat + '--:--' + tar)
        limit = len(ortholog_by_family_redundant)
        core = {}
        if not reference:
                list_ref = sorted(ortholog_by_family_redundant)
        else:
                list_ref = [reference]
                if list_ref[0] not in ortholog_by_family_redundant:
                        sys.exit('\tInvalid -reference argument: ' + ref + ' - Abort')
        print('\n#Core-genome construction...\n')
        for ref in list_ref:
                core[ref] = []
                i = len(homolog[ref][ref]) - len(ortholog_by_family_redundant[ref])
                nb_CDS = len(homolog[ref][ref])

                for CDS in ortholog_by_family_redundant[ref]:
                        i += 1
                        print('\r\t' + str(i) + '/' + str(nb_CDS) + ' CDS tested - reference: ' + ref + '\t' * 5, end = '')
                        if len(ortholog_by_family_redundant[ref][CDS]) == limit:
                                add_core = True
                                liste_w = [x.split('--:--')[0] for x in ortholog_by_family_redundant[ref][CDS]]
                                for w in ortholog_by_family_redundant[ref][CDS]:
                                        prot_w, ref_w = w.split('--:--')
                                        if ref_w in ortholog_by_family_redundant and prot_w in ortholog_by_family_redundant[ref_w] and len(ortholog_by_family_redundant[ref_w][prot_w]) == limit and (CDS + '--:--' + ref) in ortholog_by_family_redundant[ref_w][prot_w]:
                                                for e in ortholog_by_family_redundant:
                                                        if prot_w in homolog[ref][e]:
                                                                for j in ortholog_by_family_redundant[e][homolog[ref][e][prot_w]]:
                                                                        if j.split('--:--')[0] not in liste_w:
                                                                                add_core = False
                                                                                break
                                                        if not add_core:
                                                                break
                                        else:
                                                add_core = False
                                                break
                                if add_core:
                                        core[ref].append(CDS)
                print('\n')
        print('\nDone\n')
        return core


def main():
  
  args = parseArguments()
  path_blast_output = args.blast_output
  run = args.run
  path_output = args.output
  
  # if the output directory doesn't exist, create it
  os.makedirs(path_output, exist_ok = True)
  
  identity_threshold = args.identity_threshold
  alignment_threshold = args.alignment_threshold
  evalue_threshold = args.evalue_threshold
  identity_threshold = args.identity_threshold
  sequence_length_variation_threshold = args.variation_threshold
  
  ''' Homolog identification '''
  homolog = homolog_finder(path_blast_res = path_blast_output, identity_threshold = identity_threshold, alignment_threshold = alignment_threshold, evalue_threshold = evalue_threshold, sequence_length_variation_threshold = sequence_length_variation_threshold)
  with open(path_output + '/' + 'homolog_' + run + '.txt', 'w') as outfile:
    json.dump(homolog, outfile)

  ''' Core genome construction '''
  if os.path.exists(path_output + '/' + 'homolog_' + run + '.txt'):
    with open(path_output + '/' + 'homolog_' + run + '.txt', 'r') as infile:
      homolog = json.load(infile)
    core = core_finder(homolog, run, reference = '', path = path_output)
    with open(path_output + '/' + 'core_' + run + '.txt', 'w') as outfile:
      json.dump(core, outfile)
  else:
      sys.exit("No homolog json file found.")

if __name__ == "__main__":
    main()
