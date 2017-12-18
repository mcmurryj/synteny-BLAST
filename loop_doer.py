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
    """This function crawls through your faa database, keeping a rolling window of protein seqs as it goes.
    When you arrive at a protein of interest, it writes the rolling window to another faa db.
    There are known minor bugs w/r/t the edges of contigs.  Also room for speed improvements."""
    ###Format of input custom FAA DB:  >[gb_file_name] [Organism name] [Protein accession] \n protein sequence
    protein_window = deque(maxlen= feature_radius*2 +1)
    delim2         = "&"
    ###PHASE I:PARSE BLAST RESULTS; STORE THE GENOMIC INFORMATION AND BITSCORES OF HITS###
    print("Prepare to do loops........")
    protein_record = SeqIO.parse(blastdb, "fasta")
    for p in protein_record:
        protein_window.append(p)
        if len(protein_window) < (feature_radius + 1):  #This BS is to avoid corner case errors when you get a hit at beginning of a contig.
            pID = protein_window[-1].description        #Added the ".description"; don't know how this was running without error before as it would return a seqrecord.
        else:
            pID = protein_window[feature_radius].description
        if pID in ID_list:								#won't find clusters in the last radius unit of the faa file.
            the_right_contig = pID.split("|")[1]
            ID_list.remove(pID)
            for phit in protein_window:
                if len(phit.description.split("|")) == 5 :
                    species, contig, protID, annot, startstop = phit.description.split("|")
                else :
                    print("WATCH OUT!  Unsanitized fasta input with too many/too few delimiters.")
                if contig != the_right_contig:
                    pass	                           #originally did remove but deque mutation error
                else:
                    secondaryBLASThandle = open(secondary_BLAST_seq_file, "a")  #just call this before starting the loops?S
                    ungodly              = phit
                    ungodly.description  = pID + delim2 + phit.description + delim2 + "WARNING ID format is cluster_ID & cluster_member_ID & warning"
                    secondaryBLASThandle.write(">" + ungodly.description + "\n")
                    secondaryBLASThandle.write(str(ungodly.seq) + "\n")
    print("Done extracting context information!!!")
