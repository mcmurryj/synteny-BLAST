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
    	if not cluster_ID in data_box.keys():                      #swap for performance: if not cluster_ID in data_box :
    		data_box[cluster_ID]  = {}
    	for ali in blastrecord.alignments:
    		homologue_ID        = ali.hit_def.split(delim2)[1]
    		homologue_parent_ID = ali.hit_def.split(delim2)[0]
    		hitscore            = ali.hsps[0].score
    		if not cluster_member_ID in data_box[cluster_ID].keys():
    			data_box[cluster_ID][cluster_member_ID] = {}
    		data_box[cluster_ID][cluster_member_ID][homologue_ID] = {"bitscore"         : hitscore,
    		                                                         "homologue_parent" : homologue_parent_ID}
    return(data_box)
    print("Done writing data to dictionary!!!")
