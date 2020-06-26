# PARSER_blastp_muscle_prosite
This script contrasts query protein sequences to the protein sequences in GenBank files.  The aim is to find the query orthologes, create a tree and find conserved domains in the orthologue proteins. This SCRIPT does a sequential number of operations to the input:
  * Blastp: each query with all the protein sequences the genbank files
  * Muscle: generation of alignment document (.fasta) and tree file (.nw) for each query and its orthologues. Tree files can be visualized in:https://itol.embl.de/upload.cgi.
  * Domain search: uses PROSITEe database (prosite.dat file) provided by user to search domains in the genbank hits for each query.

## USAGE:
    Usage: python3 main.py -q <query_file_path> -s <subject_file_path> -c <coverage_cut_off> -i <identity_cut_off> -p <prosite_database_file_path> -o <output_path> -n <output_folder_name>

## OPTIONS:
    -q [--query]: path provided must correspond with a file with FASTA format containing query IDs and sequences. 
    -s [--subject]: path provided must be a directory containing ONLY genbank formated files. 
    -c [--qcovs]: coverage cut off must be a numeric value in 0-100. 
    -i [--pident]: identity cut off must be a numeric value in 0-100. 
    -p [--prosite_database]: path provided must correspond with prosite database file (file.dat format) 
    -o [--output]: path provided for output folder. 
    -n [--name]: output folder name will be 'parser_date_name'. the 'name' provided should not be longer than 10 characters.
All the options above must be provided

    -h [--help]: for script help information.

## INSTALLATION REQUIREMENTS:
Must have installed:
  * Biopython
  * Blast
  * Muscle
  * Pandas python module

## OUTPUT
The output is provided as a folder (folder name = 'parser_date_name'; *_\*name corresponds with the name provided by user with -n*) that contains:

  * **subject_multi.fasta** file: This is the multifasta file for allthe protein sequences for all the genbank files in the subject (-s) folder. The header format is: >protein_id@genbank_accession.
  * **organism_name_dictionary.txt** file: contains the full names of the organism of each genbank file and its genbank accession
  * **Query folders**. A folder is created per query sequence in the input query fasta file. Each query folder con contain:
  
        0_READ_ME.txt: file generated if there are no blastp hits or no domains where found with the prositer.
        1_query_name.fasta: query sequence in fasta format.
        2_blastp.tsv: blastp results for the query and the genbank multifasta.
        3_muscle_multifasta.fasta: multifasta that contains the query sequence and its orthologues (found through blastp).
        4_alingment_muscle.fasta: alinged multifasta.
        5_tree_muscle.nw: tree
        6_prosite_dom.txt: domains from query orthologues (genbank sequences) that are found in the prosite database.
