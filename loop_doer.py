#!/usr/bin/python
from Bio.Blast import NCBIXML
from collections import deque
from Bio import SeqIO

def getIDs (blastoutputfile):
    result_handle        = open(blastoutputfile, "r")
    initial_blast_record = NCBIXML.parse(result_handle)
    blastrecord          = next(initial_blast_record)
    ID_list              = []
    for ali in blastrecord.alignments:
        ID_list.append(ali.hit_def)

    return(ID_list)
    print("Parsed initial BLAST results!")

def do_loops(ID_list, blastdb, feature_radius, secondary_BLAST_seq_file):
    ###PHASE IC:  GET THE GENOME CHUNKS  (or at least the proteins) AND STORE IN FASTA FILES###
    ###Format of input custom FAA DB:  >[gb_file_name] [Organism name] [Protein accession] \n protein sequence
    protein_window = deque(maxlen= feature_radius*2 +1)
    delim2         = "&"
    ###PHASE I:PARSE BLAST RESULTS; STORE THE GENOMIC INFORMATION AND BITSCORES OF HITS###
    print("Prepare to do loops........")
    protein_record = SeqIO.parse(blastdb, "fasta")
    for p in protein_record:
    	protein_window.append(p)
    	if len(protein_window) < (feature_radius + 1):  #This BS is to avoid corner case errors when you get a hit at beginning of a contig.
    		pID = protein_window[-1]
    	else:
    		pID = protein_window[feature_radius].description
    	if pID in ID_list:								#won't find clusters in the last radius unit of the faa file.
    		the_right_contig = pID.split("|")[1]
    		ID_list.remove(pID)
    		for phit in protein_window:
    			species, contig, protID, annot, startstop = phit.description.split("|")
    			if contig != the_right_contig:
    				pass	#originally did remove but deque mutation error
    			else:
    				secondaryBLASThandle = open(secondary_BLAST_seq_file, "a")  #just call this before starting the loops?S
    				ungodly              = phit
    				ungodly.description  = pID + delim2 + phit.description + delim2 + "WARNING ID format is cluster_ID & cluster_member_ID & warning"
    				secondaryBLASThandle.write(">" + ungodly.description + "\n")
    				secondaryBLASThandle.write(str(ungodly.seq) + "\n")
    print("Done extracting context information!!!")
