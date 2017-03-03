#!/usr/bin/python
from Bio.Blast import NCBIXML
def add_annot(data_box, output) :
    delim2                                      = "&"
    cluster_ID, cluster_member_ID = "", ""
    second_BLAST_rec = NCBIXML.parse(output)
    data_box = {}
    for blastrecord in second_BLAST_rec:
    	cluster_ID        = blastrecord.query_def.split(delim2)[0]
    	cluster_member_ID = blastrecord.query_def.split(delim2)[1]
        data_box[cluster_ID][cluster_member_ID]["CDD_list"] = []
    	for ali in blastrecord.alignments:
    		CDD_ID              = ali.hit_id
    		CDD_def             = ali.hit_def
    		data_box[cluster_ID][cluster_member_ID]["CDD_list"].append(CDD_ID)
