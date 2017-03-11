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
rps_db             = var_dict["rps_db"]

###MAKE FOLDERS##
os.makedirs(output_dir)
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
import parse_BLAST2
print("Get ready to parse BLAST2!!!")
data_box = parse_BLAST2.parse(outputfile2ndBLAST)
print("Done parsing BLAST2!!!")

###PHASE IF: run rps BLAST and parse
from Bio.Blast.Applications import NcbirpsblastCommandline as rpsBLAST
#If I didn't set a db for RPS-blast, not worth trying.
if rps_db != False :
    outputfilerps = os.path.abspath(secondary_BLAST_seq_dir + "/rpsBLAST_output.XML")
    rpscmd = rpsBLAST(query  = secondary_BLAST_seq_file,
                      db     = rps_db,
                      evalue = 1E-10,
                      outfmt = 5,
                      out    = outputfilerps
                      )
    print("Running rps-BLAST")
    rpscmd()
    print("Done with RPS-blast!!!")
    from rpsBLAST_tools import add_annot
    annot_dict, annot_def_dict = add_annot(data_box, outputfilerps)
    print("Done adding domain annotations!!!")
else :
    annot_dict, annot_def_dict = {}, {}

###Write data_box to pickle
#import pickle
#pickfile = os.path.abspath(output_dir + "/data_box_pickle")
#pickle.dump(data_box, open(pickfile, "wb"))

###PHASE IIA:  PRINT TABLE WITH DATA
tab_output_dir = os.path.abspath(output_dir + "/table_output")
os.makedirs(tab_output_dir)
import print_tables as pt
print("Writing tabular output...")
pt.print_tables(data_box, tab_output_dir,
                 annot_dict     = annot_dict,
                 annot_def_dict = annot_def_dict)
import print_summary_tables as pst
pst.summarize(tab_output_dir)
print("Done writing tabular output!!!")
###PHASE IIB: MAKE NETWORK AND STORE IN XGMML
g_output_dir = os.path.abspath(output_dir + "/graph_output")
os.makedirs(g_output_dir)
import write_networks
write_networks.write_networks(data_box, g_output_dir)
