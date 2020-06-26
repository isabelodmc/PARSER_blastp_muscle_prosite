#usr/bin/env python3
# -*- coding: UTF-8 -*- 
#Isabel de Moya Clark 2020

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

######################
### SCRIPT COMENTS ###
######################

# INPUT: query sequences (fasta format), genbank files (genbank format, all 
# files in the same folder) and prosite database.
# This script does a sequencial number of operations to the input:
#	1. Blastp: each query with all the protein sequences the genbank files
#	2. Muscle: generation of alignment document (.fasta) and tree file (.nw) 
# 	for each query and its orthologues. Tree files can be visualized in:
#	https://itol.embl.de/upload.cgi. 
#	3. Domain search: uses Prosite database (prosite.dat file) provided by 
#	user to search domains	in the genbank hits for each query.

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

import sys
import getopt
import os
import re
from Bio import Seq, SeqIO

#HELP FUNCTIONS
from module0_helper import help
from module0_helper import minihelp 

import module0_helper as module0
import module1_formater as module1
import module2_blaster as module2
import module3_muscler as module3
import module4_prositer as module4



#------------------------------------------------------------------------------

#TERMINAL ARGUMENTS AND OPTIONS

args = sys.argv[1:]

#Terminal options: 
opts, args = getopt.getopt(args, "hq:s:c:i:p:o:n:", ["help", "query=", 
													"subject=", "qcovs=", 
													"pident=", 
													"prosite_database=", 
													"output=", "name="])

#Number of arguments check
if len(sys.argv)<=1: #ningún argumento proporcionado
	print("\nERROR: No arguments provided.")
	minihelp()
	sys.exit()

elif len(sys.argv)>15: #demasiados argumetnos proporcionados
	print("\nERROR: Too many arguments provided.")
	minihelp()
	sys.exit()

#Argument check and variable asingment
for opt, args in opts:

	if opt in ['-q', '--query']: #QUERY: fasta file with query sequences.
		if os.path.exists(args): #Path exists.
			if os.path.isfile(args):#File exists.
				query = args
			else:
				print("\nERROR: The path provided for the query, %s, "
					"is not a file" % args)
				minihelp()
				sys.exit()
		else:
			print("\nERROR: The path provided the query, %s, does not exist"
				"." % args)
			minihelp()
			sys.exit()
		with open(query, "r") as handle: #File has a FASTA format.
			fasta = SeqIO.parse(handle, "fasta")
			if any(fasta) == False:
				print("\nERROR: The query file provided, %s, " 
					"does not have a fasta format." % query)
				minihelp()
				sys.exit()
	
	if opt in ['-s', '--subject']: #SUBJECT: folder with genbank files.
		if os.path.exists(args): #Path exists.
			if os.path.isdir(args): #Path is a directory.
				subject = args
			else:
				print("\nERROR:The path provided for the subject, %s, "
					"does not correspond with a directory." % args)
				minihelp()
				sys.exit()
		else:
			subject = args
			print("\nERROR: The subject provided for the subject, %s, "
				"does not exist." % subject)
			minihelp()
			sys.exit()
	
	if opt in ['-c', '--qcovs']: #COVERAGE cut-off for blastp output.
		try:
			qcovs = float(args) #qcovs provided is a number. 
			if 0.0 > qcovs or qcovs > 100.0: #qcovs number is in 0-100.
				print("\nERROR: The value %s provided for qcovs is out of the"
				" range 0-100." % qcovs) 
				minihelp()
				sys.exit()
			else:
				continue
		except ValueError:		
			print("\nERROR: The value %s provided for qcovs is not a number." % args)
			minihelp()
			sys.exit()
			
	if opt in ['-i', '--pident']: #IDENTITY cut-off for blastp output.
		try:
			pident = float(args) #pident provided is a number.
			if 0.0 > pident or pident > 100.0: #pident number provided is in 0-100.
				print("\nERROR: The value %s provided for qcovs is out of the range"
					" 0-100." % pident) 
				minihelp()
				sys.exit()
			else:
				continue
		except ValueError:		
			print("\nERROR: The value %s provided for pident is not a number." % args)
			minihelp()
			sys.exit()

	if opt in ['-p', '--prosite_database']: #PROSITE DATABASE
		if os.path.exists(args): #Path provided for database exists.
			if os.path.isfile(args): #Path corresponds with a file.
				prosite_database_path = args
		else:
			print("\nERROR: The path provided for the prosite database file, %s," 
				" does not exist." % args)
			minihelp()
			sys.exit()	
		
	if opt in ['-o', '-output']: #OUTPUT FOLDER PATH
		if os.path.exists(args): #Path exists.
			if os.path.isdir(args): #Path corresponds with a directory.
				output = args
			else:
				print("\nERROR: The path provided for the output, %s, does not "
					"correspond with a directory." % args)
				minihelp()
				sys.exit()
		else:
			output = args
			print("\nERROR: The output path provided, %s, does not exist." % output)
			minihelp()
			sys.exit()

	if opt in ['-n', '--name']: #OUTPUT FOLDER NAME
		if len(args) > 10: #Name length not longer than 10 characters
			print("\nERROR: The name provided is too long. No more than 10 characters."
				" NO spaces")
			minihelp()
			sys.exit()
		else:
			name = args

	if opt in ['-h', '--help']: #HELP FUNCTION
		help()
		sys.exit()

print("\nPLEASE WAIT: PATIENCE IS THE KEY TO SUCCESS :)\n")
#------------------------------------------------------------------------------

#OUTPUT FOLDER FOR OUTPUT DOCUMENTS

from datetime import datetime
now = datetime.now()
date =  now.strftime("%d") + "_" + now.strftime("%m") + "_" + now.strftime("%Y")

name_dir_output = output + '/parser_' + date + "_" + name #Output folder name

if os.path.isdir(name_dir_output): #Checking if name exists.
	index = 0
	while True:
		index += 1 #If name exists, number is added at the end.
		new_file_name = name_dir_output + "_" + str(index)
		if os.path.isdir(new_file_name):
			continue
		else:
			name_dir_output = new_file_name
			os.makedirs(name_dir_output) #Making output file.
			break
else:
	os.makedirs(name_dir_output) #Making output file.

os.chdir(name_dir_output) #Opening output folder.

#------------------------------------------------------------------------------

#READ_ME.txt file generation.
module0.read_me(name_dir_output)

#------------------------------------------------------------------------------

#GENBANK FILES TO ONE MULTIFASTA.
#All protein sequences ad ids are extracted form all genbank files and converted
#in an individual multifasta file. File is saved in output folder.

module1.gb_to_fa(subject, name_dir_output)

#Dictionary that contains the genbank accession as key and the name of organism 
#as the code.

module1.dictionary(subject, name_dir_output)

#------------------------------------------------------------------------------

#FOLDERS FOR INDIVIDUAL QUERY
#A folder for each query in multifasta query file will be created and individual
#output for each query will be located in its own query folder.

qr=open(query, "r")
lines=qr.read().split('>')
lines = lines[1:]
lines=['>'+ seq for seq in lines]
index = 1
qr.close()	

for name in lines:
	query_name= name.split('\n')[0][1:] #Extracting sequence ID for file name
	query_dir_path = name_dir_output + "/" + str(index) + "_" + query_name
	os.makedirs(query_dir_path)
	query_file_path = query_dir_path + "/1_query_" + query_name+".fasta"
	out_file=open(query_file_path , "w")
	out_file.write(name)
	out_file.close()
	index += 1

#BLASTER-P of each query and genbank multifasta
	#blastp
	module2.blastp(query_file_path, os.path.join(
					name_dir_output, "subject_multi.fasta"), query_dir_path) 
	#blastp file filtering by IDENTITY and COVERAGE cut-offs
	module2.blastp_filtered(os.path.join(query_dir_path, "blastp_unfiltered.tsv"),
							pident, qcovs, query_dir_path)
	#Removal of unfiltered blastp file
	os.remove(os.path.join(query_dir_path, "blastp_unfiltered.tsv"))
	#Generation of the blastp file that will be saved (header is added)
	module2.blastp_to_save(os.path.join(query_dir_path, "blastp_filtered.tsv"), 
							query_dir_path)
	
	#Checking if there were any blastp hits for the query
	if os.path.getsize(os.path.join(query_dir_path, "blastp_filtered.tsv")):
	#If there are hits, a multifasta file is created for muscle analysis.
		module3.multifasta(os.path.join(query_dir_path, "blastp_filtered.tsv"), 
							query_file_path, query_dir_path)
		
#PROSITER: DOMAIN SEARCH IN PROSITE DATABASE
		module4.prosite_fasta( #Prosite fasta file qith query homologues
			os.path.join(query_dir_path, "blastp_filtered.tsv"), 
			query_file_path, query_dir_path)
		module4.prosite_parser( #Domain search
			os.path.join(query_dir_path, "prosite_fasta.fasta"), 
			prosite_database_path, query_dir_path)
			#Print message if there are no prosite hits
		if not os.path.isfile(os.path.join(query_dir_path, "6_prosite_dom.txt")):
			module4.no_prosite_hits(query_dir_path)
		#Remove files that do not need to be saved
		os.remove(os.path.join(query_dir_path, "blastp_filtered.tsv"))
		os.remove(os.path.join(query_dir_path, "prosite_fasta.fasta"))

#MUSCLER: ALINGMENT AND TREE GENERATION	
		module3.muscle_alingment(
			os.path.join(query_dir_path, "3_muscle_multifasta.fasta"),
			 query_dir_path) #ALINGMENT
		module3.muscle_maketree(
			os.path.join(query_dir_path, "4_alignment_muscle.fasta"), 
			query_dir_path) #TREE

	else: 
		#If no hits are found, the analysis will stop for the query.
		module2.no_blastp_hits(query_dir_path) #No hit file is created
		os.remove(os.path.join(query_dir_path, "2_blastp.tsv"))
		os.remove(os.path.join(query_dir_path, "blastp_filtered.tsv"))
		continue	
