#!/usr/bin/python
from Bio.Blast import NCBIXML

def parse(outputfile2ndBLAST) :
    delim2                                      = "&"
    cluster_ID, cluster_member_ID, homologue_ID = "", "", ""
    result_handle = open(outputfile2ndBLAST, "r")
    second_BLAST_rec = NCBIXML.parse(result_handle)
    data_box = {}
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
    		data_box[cluster_ID][cluster_member_ID][homologue_ID] = {"bitscore"         : hitscore,
    		                                                         "homologue_parent" : homologue_parent_ID }
            data_box[cluster_ID][cluster_member_ID]["domains"]    = []

    return(data_box)
    print("Done writing data to dictionary!!!")

def annotate(data_box, outputfileCDD) :
    delim2  = "&"
    delim = "|"

    result_handle = open(outputfileCDD, "r")
    CDD_BLAST_rec = NCBIXML.parse(outputfileCDD)

    for blastrecord in CDD_BLAST_rec:
        cluster_ID        = blastrecord.query.split(delim2)[0]
        cluster_member_ID = blastrecord.query.split(delim2)[1]
        for ali in blastrecord.alignments:
            domain_ID = ali.hit_def
            if not cluster_ID in data_box.keys():
                print("OH NO "*1000)
            else:
                data_box[cluster_ID][cluster_member_ID]["domains"].append(domain_ID)

    return(data_box)
