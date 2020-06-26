#usr/bin/env python3
# -*- coding: UTF-8 -*- 
#Isabel de Moya Clark 2020

#------------------------------------------------------------------------------
######################################################
### MUSCLER (alingment and neighbour joining tree) ###
######################################################c
#------------------------------------------------------------------------------

import sys
import os
import re
from Bio import Seq
from Bio import SeqIO
import pandas as pd

#------------------------------------------------------------------------------

def multifasta (blastp_file,  query_fasta, output_path):
	'''
	Generates multifasta file with the query sequence and its homologues
	This file will be further analyzed with muscle to obtain the alingment.
	'''
	blastp_filtered = pd.read_csv(blastp_file, header=None, 
									index_col=False, sep='\t')
	query = open(query_fasta, 'r')
	sys.stdout = open(os.path.join(output_path, 
									'3_muscle_multifasta.fasta'), 'a')
	for line in query.readlines():
		print(line.strip())
	for index, row in blastp_filtered.iterrows():
		print(">" + str(row[0]) + "\n" + str(row[4]))
	sys.stdout.close()

def muscle_alingment(multifasta, output_path):
	'''
	Generates alingment document for query and orthologues.
	'''
	command_line= ('muscle -in ' + multifasta + ' -out ' + 
					os.path.join(output_path, '4_alignment_muscle.fasta'))
	os.system(command_line)

def muscle_maketree(alignment, output_path):
	'''
	Generates neighbour joining tree document for query and orthologues.
	'''
	command_line= ('muscle -maketree -in ' + alignment + ' -out ' + 
					os.path.join(output_path, '5_tree_muscle.nw') + 
					' -cluster neighborjoining')
	os.system(command_line)
