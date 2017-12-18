#!/usr/bin/python
from Bio.Blast import NCBIXML
def add_annot(data_box, output) :
    """Parse output of RPS BLAST.
       Save the results in a dictionary.
       Output:
       annot_dict is keyed on cluster member ID, values are CDD or Pfam IDs.
       annot_def_dict is keyed on cluster member ID, values are annotations."""
    annot_dict     = {}
    annot_def_dict = {}
    rps_handle                    = open(output, "r")
    delim2                        = "&"
    cluster_ID, cluster_member_ID = "", ""
    rps_BLAST_rec = NCBIXML.parse(rps_handle)
    for blastrecord in rps_BLAST_rec:
        cluster_member_ID = blastrecord.query.split(delim2)[1]
        if not cluster_member_ID in annot_dict.keys():
            annot_dict[cluster_member_ID]    = []
            annot_def_dict[cluster_member_ID] = []
        for ali in blastrecord.alignments:
            CDD_ID              = ali.hit_id
            CDD_def             = " ".join(ali.hit_def.split(", ")[:2])
            annot_dict[cluster_member_ID].append(CDD_ID)
            annot_def_dict[cluster_member_ID].append(CDD_def)

    return(annot_dict, annot_def_dict)
