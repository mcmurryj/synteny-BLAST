#!/usr/bin/python
###INPUT USER INPUT STUFF####https://docs.python.org/2/howto/argparse.html#id1
import argparse
def get_parms():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help = "your input protein sequence in fasta format")
    parser.add_argument("database", help = "your input BLAST database")
    parser.add_argument("output", help = "the name of the directory where you wish to store output")
    args = parser.parse_args()
    return (
    "query_fa"           : args.input,
    "output_dir"         : args.output,
    "query_eval"         : 1E-30,				#cutoff evalue for the initial BLAST search
    "feature_radius"     : 5,		#how many proteins on either side to take for subsequent analysis?
    "max_no_seqs"        : 200,				#maximum number of hits to return in the initial search
    "handshake_eval"	 : 1e-20,			#cutoff evalue for the _blast_seq blast searches
    "blastdb"            : args.database,					#use nr as the db for initial blast search by default
    )
    print("Done setting inital variables.  Input sequence from: " + query_fa + "; Output goes in " + output_dir)
