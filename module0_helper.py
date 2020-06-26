#usr/bin/env python3
# -*- coding: UTF-8 -*- 
#Isabel de Moya Clark 2020

#------------------------------------------------------------------------------
#####################
### HELPER MODULE ###
#####################
#------------------------------------------------------------------------------

import sys
import os

#------------------------------------------------------------------------------

def help():
   '''
   HELP message for the script.
   '''
   print("\n###SCRIPT###")
   print("\t1. Blastp of each query from fasta file with all the protein"
      " sequences in the genbank files (2_blastp.tsv file)")
   print("\t2. MUSCLE: alingment of each query with homologues obtained from "
      "blastp (4_alingment_muscle.fasta file) and neighbor joining tree ("
      "5_tree_muscle.nw file that can be visualized in:"
      " --> https://itol.embl.de/upload.cgi).")
   print("\t3. DOMAIN SEARCH: domain search using PROSITE DATABASE in the"
      "homologues from the genbank files found through blastp.")

   print("\n###USAGE###")
   print("python3 main.py -q <query_file_path> -s <subject_file_path> "
      "-c <coverage_cut_off> -i <identity_cut_off> -p <prosite_database_file_p"
      "ath> -o <output_path> -n <output_folder_name>")

   print("\n###OPTIONS###")
   print("\t-q [--query]: path provided must correspond with a file with FASTA "
      "format containing query IDs and sequences.")
   print("\t-s [--subject]: path provided must be a directory containing ONLY"
      " genbank formated files.")
   print("\t-c [--qcovs]: coverage cut off must be a numeric value in 0-100.")
   print("\t-i [--pident]: identity cut off must be a numeric value in 0-100.")
   print("\t-p [--prosite_database]: path provided must correspond with "
      "prosite database file (file.dat format)")
   print("\t-o [--output]: path provided for output folder.")
   print("\t-n [--name]: output folder name will be 'parser_date_name'. the"
      " 'name' provided should not be longer than 10 characters.")
   print("\n##REQUIREMENTS###")
   print("Must have installed:")
   print("\t* Biopython")
   print("\t* Blast")
   print("\t* Muscle")
   print("\t* Pandas python module\n")

def minihelp():
   '''
   Simplified version of the HELP message for the script.
   Automatically printed in teh terminal if input information is wrong.
   '''
   print("\n###USAGE###")
   print("python3 main.py -q <query_file_path> -s <subject_file_path> "
      "-c <coverage_cut_off> -i <identity_cut_off> -p <prosite_database_file_p"
      "ath> -o <output_path> -n <output_folder_name>\n") 
   print("For more information: 'python3 script.py -h [--help]'.\n")

def read_me(output_path):
   '''
   Generates a file.txt that explains the contents of the output folder.
   '''
   sys.stdout = open(os.path.join(output_path, 'READ_ME.txt'), 'a')
   print("##########")
   print("##PARSER##")
   print("##########\n")
   print("Output folder contents:\n")
   print("\t* subject_multi.fasta: This is the multifasta file for all"
      "the protein sequences for all the genbank files in the subject "
      "(-s) folder. The header format is: >protein_id@genbank_accession.")
   print("\t* organism_name_dictionary.txt: contains the full names of "
      "the organism of each genbank file and its genbank accession")
   print("\t* Query folders. A folder is created per query sequence in "
      "the input query fasta file. Each query folder contains:\n")
   print("\t\t+ 0_READ_ME.txt: file generated if there are no blastp "
      "hits or no domains where found with the prositer ")
   print("\t\t+ 1_query_name.fasta: query sequence in fasta format")
   print("\t\t+ 2_blastp.tsv: blastp results for the query and the "
      "genbank multifasta")
   print("\t\t+ 3_muscle_multifasta.fasta: multifasta that contains "
      "the query sequence and its orthologues (found through blastp).")
   print("\t\t+ 4_alingment_muscle.fasta: alinged multifasta.")
   print("\t\t+ 5_tree_muscle.nw: tree")
   print("\t\t+ 6_prosite_dom.txt: domains from query orthologues "
      "(genbank sequences) that are found in the prosite database.")
   sys.stdout.close()
