#usr/bin/env python3
# -*- coding: UTF-8 -*- 
#Isabel de Moya Clark 2020

#------------------------------------------------------------------------------
########################
### BLASTER (blastp) ###
########################
#-------------------------------------------------------------------------------

import sys
import os
import csv
import numpy as np
import pandas as pd

def blastp(query, subject, output_path):
	'''
	Does blastp with a given query and subject. 
	Output file path must be provided.
	'''
	command_line = ('blastp -query ' + query + ' -subject ' + subject + 
	' -evalue 0.00001 -outfmt "6 sseqid pident qcovs evalue sseq" > ' +
	os.path.join(output_path, 'blastp_unfiltered.tsv'))

	os.system(command_line)


def blastp_filtered(blastp_file, identity, coverage, output_path):
	'''
	Generates a blastp file that is filtered by identity and coverage cutoffs 
	and that contains a header.
	'''
	blastp_unfiltered = pd.read_csv(blastp_file, header=None, 
									index_col=False, sep='\t')
	sys.stdout = open(os.path.join(output_path, 'blastp_filtered.tsv'), 'a')
	
	for index, row in blastp_unfiltered.iterrows():
	
		if float(row[1]) >= identity and float(row[2]) >= coverage:
			print(str(row[0]) + "\t" + str(row[1]) + 
					"\t" + str(row[2]) + "\t" + str(row[3]) + "\t" + str(row[4]))

def blastp_to_save (blastp_file, output_path):
	'''
	Generates a blastp file with a header that will be saved
	in the output file query folder.
	'''
	blastp_filtered = open(blastp_file, 'r')
	sys.stdout = open(os.path.join(output_path, '2_blastp.tsv'), 'a')
	print("sseqid\tpident\tqcovs\tevalue\tsseq")
	
	for line in blastp_filtered.readlines():
		print(line.strip())

def no_blastp_hits(output_path):
	'''
	If no blastp hits are found for a query, the program for that 
	query will not continue and a file with "NO HITS" message will 
	be created.
	'''
	sys.stdout = open(os.path.join(output_path, '0_READ_ME.txt'), 'a')
	print("NO bastp HITS for this query.")
