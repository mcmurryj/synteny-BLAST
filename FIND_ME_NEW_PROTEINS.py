###PHASE -1 : IMPORT LIBRARIES###
#!/usr/bin/python
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
import os

###INPUT USER INPUT STUFF####https://docs.python.org/2/howto/argparse.html#id1
import get_parms
var_dict           = get_parms.get_parms()
query_fa           = var_dict["query_fa"]
output_dir         = os.path.abspath(var_dict["output_dir"])
query_eval         = var_dict["query_eval"]				#cutoff evalue for the initial BLAST search
feature_radius     = var_dict["feature_radius"]			#how many proteins on either side to take for subsequent analysis?
handshake_eval	   = var_dict["handshake_eval"]			#cutoff evalue for the _blast_seq blast searches
blastdb            = var_dict["blastdb"]					#use nr as the db for initial blast search by default

###MAKE FOLDERS##
from make_folders import make_folders
make_folders(output_dir)
blastoutputdir    = os.path.abspath(output_dir + "/initial_BLAST_data")
os.makedirs(blastoutputdir)
blastoutputfile   = os.path.abspath(blastoutputdir + "/initial_blast_output.XML")

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
os.makedirs(secondary_BLAST_seq_dir)
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
<<<<<<< HEAD
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
###PHASE IIA1:  Run HMMscan on the BLAST db and parse results
print("Running HMMscan on the secondary blast DB.....")
from hmmscan_command import run_hmmscan
run_hmmscan(fasta = secondary_BLAST_seq_file)
print("Parsing hmmscan results.......")
from hmmscan_parse import hmmscan_parse
pfam_dict                = hmmscan_parse(os.path.abspath("./hmmscan_data/hmmscan_out.tab"))   #hardcoded path is aesthetically unappealing
cluster_member_pfam_dict = pfam_dict["cluster member ID pfam dict"]
print("Done recording dict of pfam domains!!!")

###PHASE IIA:  PRINT TABLE WITH DATA
delim              = "|"
header             = ["cluster member ID", "cluster member annotation", "cluster member startstop",
                      "homologue organism", "homologue contig", "homologue proteinID", "homologue annot",
					  "homologue startstop", "BLAST bitscore", "cluster member pfam domains"]
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
				data2     = homologue_ID.split(delim)
				score     = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
				pfam_hits = cluster_member_pfam_dict[cluster_member_ID]
				values = data1[2:] + data2 + score + pfam_hits				#will make ragged table
				tabout_handle.write('\t'.join(map(str, values))+"\n")
			else :                                       #less interesting-where a protein is related to itself.
				score = [data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]]
				values = data1[2:] + ["NA"]*5 + score + pfam_hits
				tabout_handle.write('\t'.join(map(str, values))+"\n")
=======
import parse_BLAST2
data_box = parse_BLAST2.parse(outputfile2ndBLAST)

###PHASE IIA:  PRINT TABLE WITH DATA
tab_output_dir = os.path.abspath(output_dir + "/table_output")
os.makedirs(tab_output_dir)
import print_tables
print_tables.print_tables(data_box, tab_output_dir)
>>>>>>> ae02f17b845375a777a12359230557d0f6304f7e

###PHASE IIB: MAKE NETWORK AND STORE IN XGMML
g_output_dir = os.path.abspath(output_dir + "/graph_output")
os.makedirs(g_output_dir)
import write_networks
write_networks.write_networks(data_box, g_output_dir)
