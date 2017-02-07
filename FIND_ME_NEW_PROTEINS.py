#!/usr/bin/python
###PHASE -1 : IMPORT LIBRARIES###
from sys import argv
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
import shutil
import glob
from collections import deque
import operator

###PHASE 0:  PRINT OUT A HELPFUL HELP STATEMENT AND GET INPUTS##
###INPUT USER INPUT STUFF####https://docs.python.org/2/howto/argparse.html#id1
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input", help = "your input protein sequence in fasta format")
parser.add_argument("database", help = "your input BLAST database")
parser.add_argument("output", help = "the name of the directory where you wish to store output")
args = parser.parse_args()
query_fa           = args.input
output_dir         = args.output
query_eval         = 1E-30				#cutoff evalue for the initial BLAST search
feature_radius     = 5			#how many proteins on either side to take for subsequent analysis?
max_no_seqs        = 200				#maximum number of hits to return in the initial search
handshake_eval	   = 1e-20			#cutoff evalue for the secondary blast searches
blastdb            = args.database					#use nr as the db for initial blast search by default
print("Done setting inital variables.  Input sequence from: " + query_fa + "; Output goes in" + output_dir)

###MAKE DATA STRUCTURES###sillycomment
data_box = {}	#databox = {clusterID : {clustermemberID : {homologueID : bitscore, .....}, ......}, ......}

###MAKE FOLDER FOR OUTPUT###
if os.path.exists(output_dir):
	shutil.rmtree(output_dir)
os.makedirs(output_dir)

###MAKE SUBFOLDER FOR BLAST OUTPUT####
blastoutputdir = output_dir + "/"+ "initial_BLAST_data"
os.makedirs(os.path.abspath(blastoutputdir))
blastoutputfile = os.path.abspath(blastoutputdir + "/" + "initial_blast_output.XML")

###PHASE IA: BLAST QUERY SEQUENCE AGAINST DATABASE###
cline = NcbiblastpCommandline(query = query_fa, db = blastdb, evalue = query_eval, outfmt = 5, out = blastoutputfile)
cline()				#to actually run the thing
print("Completed initial BLAST search")

###PHASE I:PARSE BLAST RESULTS; STORE THE GENOMIC INFORMATION AND BITSCORES OF HITS###
result_handle = open(blastoutputfile, "r")
initial_blast_record = NCBIXML.parse(result_handle)
blastrecord = next(initial_blast_record)			  					#are they ordered best hit first?  THis is needed.  lAso that the query is in refseq NR.
ali               = blastrecord.alignments[0]											#best hit first
query_ID          = ali.hit_def									#get the ID of the presumptive query  Will need work on ID retrieval.
cluster_ID        = query_ID
cluster_member_ID = query_ID
homologue_ID      = query_ID
result_handle.close()
print("Parsed initial BLAST results!")

data_box[cluster_ID] = {cluster_member_ID : {}}
####THis chunk has some issues.
for ali in blastrecord.alignments:
	homologue_ID = ali.hit_def	 													# I forget what alignments vs. hsps are.  XXX.title or .hit_def is more info, .accession is like the WP ####
	hitscore = ali.hsps[0].score													#I think this is OK?
	data_box[cluster_ID][cluster_member_ID][homologue_ID] = hitscore			#strictly speaking program will just write over this later.
print("Read in scores of best hits to query sequence!")

###PHASE IC:  GET THE GENOME CHUNKS  (or at least the proteins) AND STORE IN FASTA FILES###
###Format of input custom FAA DB:  >[gb_file_name] [Organism name] [Protein accession] \n protein sequence
###make some data structures
protein_window = deque(maxlen= feature_radius*2 +1)
temp_ID_list = list(data_box[query_ID][query_ID].keys())
delim2 = "&"

###make output folder
secondary_BLAST_seq_dir = output_dir + "/" + "BLAST2_data"
os.makedirs(os.path.abspath(secondary_BLAST_seq_dir))
secondary_BLAST_seq_file = os.path.abspath(secondary_BLAST_seq_dir + "/" + "secondary_BLAST_seqs.faa")
print("Made dir " + secondary_BLAST_seq_dir + " to hold file " + secondary_BLAST_seq_file)

###do loops
print("Prepare to do loops........")
protein_record = SeqIO.parse(blastdb, "fasta")
for p in protein_record:
	protein_window.append(p)								###again I will need to  harmonize the ID system
	if len(protein_window) < (feature_radius + 1):
		pID = protein_window[-1]
	else:
		pID = protein_window[feature_radius].description
	if pID in temp_ID_list:								#won't find clusters in the last radius unit of the faa file.
		the_right_contig = pID.split("|")[1]
		temp_ID_list.remove(pID)
		for phit in protein_window:
			species, contig, protID, annot, startstop = phit.description.split("|")
			if contig != the_right_contig:
				pass	#originally did remove but deque mutation error
			else:
				secondaryBLASThandle = open(secondary_BLAST_seq_file, "a")
				ungodly = phit
				ungodly.description = pID + delim2 + phit.description + delim2 + "WARNING ID format is cluster_ID & cluster_member_ID & warning"
				SeqIO.write(ungodly, secondaryBLASThandle, 'fasta')
print("Done extracting context information!!!")

###ALSO WRITE NEW BLAST DB###
make2ndBLASTdbcmd = "makeblastdb -in " + secondary_BLAST_seq_file + " -input_type fasta -dbtype prot"
print(make2ndBLASTdbcmd)
subprocess.call(make2ndBLASTdbcmd, shell = True)
print("Done making secondary BLAST DB!!!")

###PHASE ID:  run 2ndary blast
outputfile2ndBLAST = os.path.abspath(secondary_BLAST_seq_dir + "/" + "2ndBLAST_output.XML")
run2ndBLASTdbcmd = NcbiblastpCommandline(
query = secondary_BLAST_seq_file,
db = secondary_BLAST_seq_file,
evalue = handshake_eval,
outfmt = 5,
out = outputfile2ndBLAST)
print("RUnning command: " + str(run2ndBLASTdbcmd))
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
delim = "|"
header = ["organism", "contig", "cluster member ID", "cluster member annotation", "cluster member startstop",
"homologue organism", "homologue contig", "homologue proteinID", "homologue annot", "homologue startstop", "BLAST bitscore"]
cluster_ID, cluster_member_ID, homologue_ID = "", "", ""

tab_output_dir = os.path.abspath(output_dir + "/"+ "table_output")
os.makedirs(tab_output_dir)

for cluster_ID in data_box.keys():
	data0 = cluster_ID.split(delim)
	WP_no = data0[2]
	tabout_handle = open(os.path.abspath(tab_output_dir + "/"+ WP_no + ".tsv"), "a")
	tabout_handle.write('\t'.join(map(str,header))+"\n")
	for cluster_member_ID in data_box[cluster_ID].keys():
		data1 = cluster_member_ID.split(delim)
		for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
			if cluster_member_ID != homologue_ID:
				data2 = homologue_ID.split(delim)
				score = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
				values = data0[0:2] + data1[2:] + data2 + score
				tabout_handle.write('\t'.join(map(str, values))+"\n")
			else :
				score = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
				values = data0[0:2] + data1[2:] + ["NA"]*5 + score
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
