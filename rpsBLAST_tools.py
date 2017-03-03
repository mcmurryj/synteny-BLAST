#!/usr/bin/python
from Bio.Blast import NCBIXML
def add_annot(data_box, output) :
    annot_dict = {}
    rps_handle                    = open(output, "r")
    delim2                        = "&"
    cluster_ID, cluster_member_ID = "", ""
    rps_BLAST_rec = NCBIXML.parse(rps_handle)
    for blastrecord in rps_BLAST_rec:
        cluster_member_ID = blastrecord.query.split(delim2)[1]
        if not cluster_member_ID in annot_dict.keys():
            annot_dict[cluster_member_ID] = []
        for ali in blastrecord.alignments:
            CDD_ID              = ali.hit_id
            CDD_def             = ali.hit_def
            annot_dict[cluster_member_ID].append(CDD_ID)

    return(annot_dict)
