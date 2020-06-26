#usr/bin/env python3
# -*- coding: UTF-8 -*- 
#Isabel de Moya Clark 2020

#------------------------------------------------------------------------------
########################################
### FORMATER (genbank to multifasta) ###
########################################
#------------------------------------------------------------------------------

import sys
import os
from Bio import Seq, SeqIO
import pickle

#------------------------------------------------------------------------------

def gb_to_fa(input_path, output_path):
	'''
	This function gets all protein sequences from genbank files in a folder
	and converts them in an individual MULTIFASTA file that can be used in
	a blastp search.
	Each sequence is acmopanied by its protein ID and genome record name.
	'''	
	with open(os.path.join(output_path,"subject_multi.fasta"), "a") as out_handle:
	
		for file in os.listdir(input_path):
	
				with open(os.path.join(input_path, file), "r") as in_handle:

					for record in SeqIO.parse(in_handle, "genbank"):

							for feature in record.features:
								
								if feature.type == 'CDS':

									try:
										out_handle.write(">%s@%s\n%s\n" % (
											feature.qualifiers['locus_tag'][0], 
											record.name, 
											feature.qualifiers['translation'][0])
											)
									except:
										pass

def dictionary(genbank_directory, output_path):
	'''
	Generates a dictionary containing genbank record names as keys
	and organism names as codes for all genbank files in genbank directory
	'''
	dic_id_org={}
	
	for file in os.listdir(genbank_directory):
	
		with open(os.path.join(genbank_directory, file), "r") as in_handle:
	
			for record in SeqIO.parse(in_handle, "genbank"):
				dic_id_org[record.name] = record.annotations['organism']
	sys.stdout = open(os.path.join(output_path, "organism_name_dictionary.txt"), "w")
	print(dic_id_org)
	sys.stdout.close()

#Estuve intentando usar el diccionario para cambiar el genbank accession al 
#nombre del organismo en los archivos multifasta usados para generar el alingment 
#y el árbol pero no lo consegí así que solo he generado un archivo que contiene 
#el diccionario en formato texto y lo de abajo es el código que inteté usar para 
#que en los archivos multifasta cambiara el accession por el nombre del organismo.

#	dic_file=open(os.path.join(output_path, "dictionary.pkl"), "wb")
#	pickle.dump(dic_id_org, dic_file)
#	dic_file.close()

#def id_org_converter(dictionary_file, muscle_multifasta_file, output_path):
#	dic_file=open(dictionary_file, "rb")
#	dic=pickle.load(dic_file)
#	dic_file.close()
#
#	with open(muscle_multifasta_file, "r") as f:
#		for line in f.read().split('\n'):
#			if line.startswith(">"):
#				try:
#					prot_id = line.split("@")[0]
#					acession = line.split("@")[1:]
#					org = dic[str(accession)]
#					print (str(prot_id) + "@" + str(org))
#				except:
#					print(line)
#			else:
#				print(line)
#		else:
#			print(line)