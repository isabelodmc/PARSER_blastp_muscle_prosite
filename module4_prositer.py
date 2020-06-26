#usr/bin/env python3
# -*- coding: UTF-8 -*- 
#Isabel de Moya Clark 2020

#------------------------------------------------------------------------------
################
### PROSITER ###
################
#------------------------------------------------------------------------------

import sys
import os
import re
from Bio import Seq, SeqIO
from Bio.ExPASy import Prosite,Prodoc
import pandas as pd

#------------------------------------------------------------------------------

def prosite_fasta(blastp_file,  query_fasta, output_path):
	'''
	Generates a fasta file that contains only query orthologues.
	'''
	blastp_filtered = pd.read_csv(blastp_file, header=None, 
									index_col=False, sep='\t')
	sys.stdout = open(os.path.join(output_path, 'prosite_fasta.fasta'), 'a')
	
	for index, row in blastp_filtered.iterrows():
		print(">" + str(row[0]) + "\n" + str(row[4]))

def prosite_parser(query_path, prosite_domains, output_path):
	'''
	Finds conserved domains of query sequences (genbank orthologues) in the 
	prosite database.
	'''
	with open(query_path, "r") as qr_handle: #opening file with orthologues

		for qr_record in SeqIO.parse(qr_handle, "fasta"): 
			qr_fasta = str(qr_record.seq) 

			with open(prosite_domains,"r") as pro_handle: #opening database
				pro_records = Prosite.parse(pro_handle)

				for pro_record in pro_records:
					#establishing pattern for eficient domain search
					prot_dom = pro_record.pattern.replace("-","").replace("{","[^").replace("}","]").replace("(","{").replace(")","}").replace(".","").replace("x",".")
					prot_doms = str(prot_dom)

					if prot_doms != "" and re.search(prot_doms, qr_fasta):
						#saving relevant information for the domains found.
						with open(os.path.join(output_path, "6_prosite_dom.txt"),
								 "a") as domain_info:
							domain_info.write("\nQUERY ORTHOLOGUE: "
												+qr_record.id)
							domain_info.write("\n\t* Name: " 
												+ pro_record.name)
							domain_info.write("\n\t* Accession: " 
												+ pro_record.accession)
							domain_info.write("\n\t* Description: " 
												+ pro_record.description)
							domain_info.write("\n\t* Pattern: " 
												+ prot_doms)
							domain_info.write("\n\t* Domain: " 
												+ str(re.findall(prot_doms, qr_fasta)[0]) + "\n")

def no_prosite_hits(output_path):
	'''
	If no blastp hits are found for a query, the program for that 
	query will not continue and a file with "NO HITS" message will 
	be created.
	'''
	sys.stdout = open(os.path.join(output_path, '0_READ_ME.txt'), 'a')
	print("NO prsotie domain HITS for the homologues of this query.")
