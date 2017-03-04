#!/usr/bin/python
###INPUT USER INPUT STUFF####https://docs.python.org/2/howto/argparse.html#id1
import argparse
def get_parms():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input",
                        help = "your input protein sequence in fasta format")
    parser.add_argument("-database",
                        help = "your input BLAST database")
    parser.add_argument("-output",
                        help = "the name of the directory where you wish to store output")
    parser.add_argument("--query_eval",
                       default = 1E-30,
                        help = "The e value cutoff to use for the initial blast search")
    parser.add_argument("--feature_radius",
                         default = 5,
                         help = "How wide a window around your initial hit do you want to consider, in coding sequences?")
    parser.add_argument("--handshake_eval",
                        default = 1E-20,
                        help = "E value cutoff for second round of BLAST.")
    parser.add_argument("--rps_db",
                        default = False,
                        help = "Database of domains to use for RPS-blast, null by default")
    parser.add_argument("--random",
                        default = False,
                        help = "Accept probability for random sampling of 1st BLAST results.")

    args = parser.parse_args()
    return ({
    "query_fa"           : args.input,
    "output_dir"         : args.output,
    "query_eval"         : args.query_eval,				#cutoff evalue for the initial BLAST search
    "feature_radius"     : args.feature_radius,		 #how many proteins on either side to take for subsequent analysis?
    "max_no_seqs"        : 200,				#maximum number of hits to return in the initial search
    "handshake_eval"	 : 1e-20,			#cutoff evalue for the _blast_seq blast searches
    "blastdb"            : args.database,
    "rps_db"             : args.rps_db,
    "random"             : args.random
    })
