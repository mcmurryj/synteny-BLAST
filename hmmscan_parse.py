#!/usr/bin/python
def hmmscan_parse (hmmscan_output):
    """This function takes as input the path to an hmmscan output file.
       It returns a dict linking cluster_member_ID to a list of pfam domains
       and a dict linking the cluster ID to a list of pfam domains."""

    hmmscan_file   = os.path.abspath(hmmscan_output)
    hmmscan_handle = open(hmmscan_file, "r")
    from Bio.SearchIO import HmmerIO.hmmer3_tab.Hmmer3TabParser
    hmmscan_output = HmmerIO.hmmer3_tab.Hmmer3TabParser(HMMhandle)

    cluster_ID            = ""
    cluster_member_ID     = ""
    homologue_ID          = ""

    clustermemberID_pf_dict = {}
    clusterID_pf_dict  = {}

    for hso in hmmscan_output:
        IDlist                         = hso.QueryResult.id.split(delim2)
        cluster_member_ID              = IDlist[1]
        pfID                           = hso.hit.id
        if not cluster_member_ID in clustermemberID_pf_dict.keys() :
            clustermember_pf_dict[cluster_member_ID] = []
        clustermember_pf_dict[cluster_member_ID].append(pfID)
        if not cluster_ID in clusterID_pf_dict.keys() :
            clusterID_pf_dict[cluster_ID] = []
        clusterID_pf_dict[cluster_ID].append(pfID)

    return({"cluster member ID pfam dict" : clustermemberID_pf_dict, "cluster ID pfam dict" : clusterID_pf_dict)
