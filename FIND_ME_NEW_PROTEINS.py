#!/usr/bin/python
###PHASE -1 : IMPORT LIBRARIES###
from sys import argv
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import shutil
import glob
import operator
import get_parms()

###INPUT USER INPUT STUFF####https://docs.python.org/2/howto/argparse.html#id1
var_dict           = get_parms
query_fa           = var_dict["query_fa"]
output_dir         = os.path.abspath(var_dict["output_dir"])
query_eval         = var_dict["query_eval"]				#cutoff evalue for the initial BLAST search
feature_radius     = var_dict["feature_radius"]			#how many proteins on either side to take for subsequent analysis?
handshake_eval	   = var_dict["handshake_eval"]			#cutoff evalue for the _blast_seq blast searches
blastdb            = var_dict["blastdb"]					#use nr as the db for initial blast search by default

###MAKE FOLDERS##
from make_folders import make_folders
make_folders(output_dir)
blastoutputdir    = path.abspath(output_dir + "/initial_BLAST_data")
make_folders(blastoutputdir)
blastoutputfile = path.abspath(blastoutputdir + "/initial_blast_output.XML")

###PHASE IA: BLAST QUERY SEQUENCE AGAINST DATABASE###
cline = NcbiblastpCommandline(query  = query_fa,
                              db     = blastdb,
							  evalue = query_eval,
							  outfmt = 5,
							  out    = blastoutputfile)
cline()														#to actually run the thing
print("Completed initial BLAST search")

###make output folder for BLAST2
secondary_BLAST_seq_dir = output_dir + "/BLAST2_data"
make_folders(secondary_BLAST_seq_dir)
secondary_BLAST_seq_file = os.path.abspath(secondary_BLAST_seq_dir + "/secondary_BLAST_seqs.faa")
print("Made dir " + secondary_BLAST_seq_dir + " to hold file " + secondary_BLAST_seq_file)

###RETRIEVE
import loop_doer
id_list  = loop_doer.getIDs(blastoutputfile)
loop_doer.do_loops(id_list, blastdb, feature_radius, secondary_BLAST_seq_file)

###ALSO WRITE NEW BLAST DB###
make2ndBLASTdbcmd = "makeblastdb -in " + secondary_BLAST_seq_file + " -input_type fasta -dbtype prot"
print(make2ndBLASTdbcmd)
subprocess.call(make2ndBLASTdbcmd, shell = True)
print("Done making secondary BLAST DB!!!")

###PHASE ID:  run 2ndary blast
outputfile2ndBLAST = os.path.abspath(secondary_BLAST_seq_dir + "/2ndBLAST_output.XML")
run2ndBLASTdbcmd = NcbiblastpCommandline(
query   = secondary_BLAST_seq_file,
db      = secondary_BLAST_seq_file,
evalue  = handshake_eval,
outfmt  = 5,
out     = outputfile2ndBLAST)
print("Running command: " + str(run2ndBLASTdbcmd))
run2ndBLASTdbcmd()
print("Done running secondary BLAST search!!!")

###Phase 1E:PARSE BLAST OUTPUT; STORE BITSCORES IN data_box[initialhitID][clustermemberID][homologueID] = bitscore
###AM HAVING PROBLEMS STORING DATA IN DICT WITHOUT OVERWRITING PREVIOUS ENTRIES
cluster_ID, cluster_member_ID, homologue_ID = "", "", ""
result_handle = open(outputfile2ndBLAST, "r")
second_BLAST_rec = NCBIXML.parse(result_handle)
for blastrecord in second_BLAST_rec:
	cluster_ID        = blastrecord.query.split(delim2)[0]
	cluster_member_ID = blastrecord.query.split(delim2)[1]
	if not cluster_ID in data_box.keys():
		data_box[cluster_ID]  = {}
	for ali in blastrecord.alignments:
		homologue_ID        = ali.hit_def.split(delim2)[1]
		homologue_parent_ID = ali.hit_def.split(delim2)[0]
		hitscore            = ali.hsps[0].score
		if not cluster_member_ID in data_box[cluster_ID].keys():
			data_box[cluster_ID][cluster_member_ID] = {}
		data_box[cluster_ID][cluster_member_ID][homologue_ID] = {"bitscore" : hitscore,
		                                                         "homologue_parent" : homologue_parent_ID}
print("Done writing data to dictionary!!!")

###PHASE IIA:  PRINT TABLE WITH DATA
delim              = "|"
header             = ["cluster member ID", "cluster member annotation", "cluster member startstop",
"homologue organism", "homologue contig", "homologue proteinID", "homologue annot", "homologue startstop", "BLAST bitscore"]
cluster_ID         = ""
cluster_member_ID  = ""
homologue_ID       = ""

#Make the directory
tab_output_dir = os.path.abspath(output_dir + "/table_output")
os.makedirs(tab_output_dir)

#iterate thru all the info and print stuff out
for cluster_ID in data_box.keys():
	#get data chunks from the titles of the fasta sequences
	data0         = cluster_ID.split(delim)
	WP_no         = data0[2]				#third element is the WP_number
	org_cont      = data0[0:2]				#first 2 elements are species and contig
	tabout_handle = open(os.path.abspath(tab_output_dir + "/"+ WP_no + ".tsv"), "a")
	tabout_handle.write('\t'.join(map(str, org_cont)) + "\n")		#print species and contig at the top; these are invariant
	tabout_handle.write('\t'.join(map(str, header  )) + "\n")		#print the header info
	for cluster_member_ID in data_box[cluster_ID].keys():
		data1 = cluster_member_ID.split(delim)
		for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
			if cluster_member_ID != homologue_ID:       #if it is not self-identity
				data2 = homologue_ID.split(delim)
				score = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
				values = data1[2:] + data2 + score
				tabout_handle.write('\t'.join(map(str, values))+"\n")
			else :                                       #less interesting-where a protein is related to itself.
				score = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
				values = data1[2:] + ["NA"]*5 + score
				tabout_handle.write('\t'.join(map(str, values))+"\n")

###PHASE IIB: MAKE NETWORK AND STORE IN XGMML
import networkx as nx
cluNet = nx.Graph()
protNet = nx.Graph()
cluster_ID, cluster_member_ID, homologue_ID, homologue_parent_ID = "", "", "", ""

for cluster_ID in data_box.keys():
	cluNet.add_node(cluster_ID)

for cluster_ID in data_box.keys():
	for cluster_member_ID in data_box[cluster_ID].keys():
		for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
			score               = data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]
			homologue_parent_ID = data_box[cluster_ID][cluster_member_ID][homologue_ID]["homologue_parent"]
			IDtup               = tuple(sorted([cluster_ID, homologue_parent_ID]))
			IDtup2              = tuple(sorted([cluster_member_ID, homologue_ID]))
			protNet.add_node(IDtup2[0])
			protNet.add_node(IDtup2[1])
			protNet.add_edge(IDtup2[0], IDtup2[1], weight = score)
			if IDtup in cluNet.edges():
				cluNet[IDtup[0]][IDtup[1]]["weight"] += score
			else:
				cluNet.add_edge(IDtup[0], IDtup[1], weight = score)

#tidy up the network by removing self-connections and deleting loner nodes.
for n in protNet.nodes():
	if len(protNet[n].keys()) <= 1:
		protNet.remove_node(n)
	if protNet.has_edge(n, n):
		protNet.remove_edge(n,n)

for n in cluNet.nodes():
	if cluNet.has_edge(n, n):
		cluNet.remove_edge(n,n)

g_output_dir = os.path.abspath(output_dir + "/"+ "graph_output")
os.makedirs(g_output_dir)
nx.write_gml(cluNet, os.path.abspath(g_output_dir + "/"+ "graph_o_clusters.gml"))
nx.write_gml(protNet, os.path.abspath(g_output_dir + "/"+ "graph_o_proteins.gml"))

		#black magic-returns list of tuples sorted by the value of the second element of the tuple
		#in this case that element is the bitscore of the hit.
		#Print the SECOND best one.  Best one is self.
		#entry_list = sorted(data_box[cluster_ID][cluster_member_ID].items(), key=operator.itemgetter(1))
		#if len(entry_list) == 1:
		#	data2 = ["NO HITS BESIDES SELF"]
		#	score = ["NA"]
		#else:
		#	data2 =entry_list[1][0].split(delim)
		#	score = [str(entry_list[1][1])]
		#print(','.join(map(str, header)))
		#print(','.join(map(str, data0 + data1 + data2 + score)))
###PHASE IIB: MAKE EDGE-WEIGHTED NETWORK
###PHASE IIC###
###ANNOTATE EACH PROTEIN WITH PFAM/COG ASSIGNMENTS###
###FOR EACH INITIAL HIT, PRINT TABLE WITH COLUMNS:  PROTEIN ID/PFAM/COG/EVAL
#Species, Contig, Accession, Start-Stop (BP), Annotation, Cluster ID, Best hit ID, Best hit score
