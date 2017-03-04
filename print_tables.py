#!/usr/bin/python
#TODO:  change the tabular output so it is sorted on homologue_ID then cluster_member_ID
import os
def print_tables(data_box, tab_output_dir, annot_dict = {}, annot_def_dict = {}):
	"""Iterate thru all the info in data_box and print to tab_output_dir"""

	header = ["cluster member ID", "cluster member annotation", "cluster member startstop",
		      "homologue organism", "homologue contig", "homologue proteinID", "homologue annot",
			  "homologue startstop", "BLAST bitscore", "CDD domains", "Domain definitions"]
	delim  = "|"

	for cluster_ID in data_box.keys():
		#get data chunks from the titles of the fasta sequences
		data0         = cluster_ID.split(delim)
		WP_no         = data0[2]				#third element is the WP_number
		org_cont      = data0[0:2]				#first and second elements are the species, contig
		#Open a file handle and write header information
		tabout_handle = open(os.path.abspath(tab_output_dir + "/"+ WP_no + ".tsv"), "a")
		tabout_handle.write('\t'.join(map(str, org_cont)) + "\n")		#print species and contig at the top; these are invariant
		tabout_handle.write('\t'.join(map(str, header  )) + "\n")		#print the header info
		#Loop thru the dict and print out the good stuff.
		for cluster_member_ID in data_box[cluster_ID].keys():
			data1 = cluster_member_ID.split(delim)
			#Check to see if there is a domain annotation list; if there is add it.
			if cluster_member_ID in annot_dict.keys() :
				annot_list     = annot_dict[cluster_member_ID]
				annot_def_list = annot_def_dict[cluster_member_ID]
			else :
				annot_list     = []
				annot_def_list = []
			for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
				if cluster_member_ID != homologue_ID :       #if it is not self-identity
					data2 = homologue_ID.split(delim)
					score = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
					annot = [",".join(annot_list), ",".join(annot_def_list)]
					values = data1[2:] + data2 + score + annot
					tabout_handle.write('\t'.join(map(str, values))+"\n")
				else :                                       #less interesting-where a protein is related to itself.
					score  = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
					annot  = [",".join(annot_list), ",".join(annot_def_list)]
					values = data1[2:] + ["NA"]*5 + score + annot
					tabout_handle.write('\t'.join(map(str, values))+"\n")
