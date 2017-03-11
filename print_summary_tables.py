import pandas as pd
from os.path import abspath
from glob import glob
def summarize(tab_output_dir) :
    filelist  = glob(abspath(tab_output_dir + "/*.tsv"))
    for f in filelist :
        tabdata   = pd.read_csv(f, header = 1, sep = "\t")
        adict = {}
        for i in set(tabdata["cluster member ID"]) :
            goodrows    = tabdata[tabdata["cluster member ID"] == i]
            hitcounts   = int(len(goodrows["cluster member ID"]))
            cdd_doms    = goodrows["CDD domains"].iloc[0]
            cdd_counts  = int(sum(tabdata["CDD domains"] == cdd_doms))
            cdd_des     = goodrows["Domain definitions"].iloc[0]
            adict[i]    = {"No._of_homologues"       : hitcounts,
                           "No._of_same_pfam_domain" : cdd_counts,
                           "Pfam_domain_def."        : cdd_des}
            tabsum      = pd.DataFrame(adict).transpose()
            tabsum.sort_values('No._of_homologues', ascending = False, inplace = True)
            tabsum.index.rename("NCBI protein ID", inplace = True)
            tabsum.to_csv(abspath(f[:-4] + "synop.tsv"), sep = "\t")

def print_tables(data_box, tab_output_dir, annot_dict = {}, annot_def_dict = {}):
    """Iterate thru all the info in data_box and print to tab_output_dir"""
    delim  = "|"
    for cluster_ID in data_box.keys():
        cdict         = {}
        #get data chunks from the titles of the fasta sequences
        data0         = cluster_ID.split(delim)
        species       = data0[0].replace(" ", "_")
        WP_no         = data0[2]                #third element is the WP_number
        org_cont      = data0[0:2]                #first and second elements are the species, contig
        idx           = 0
        #Loop thru the dict and print out the good stuff.
        for cluster_member_ID in data_box[cluster_ID].keys():
            data1 = cluster_member_ID.split(delim)
            #Check to see if there is a domain annotation list; if there is add it.
            if cluster_member_ID in annot_dict.keys() :
                annot_list     = annot_dict[cluster_member_ID]
                annot_def_list = annot_def_dict[cluster_member_ID]
            else :
                annot_list     = []
                annot_def_list = []
            for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
                if cluster_member_ID != homologue_ID :       #if it is not self-identity
                    data2 = homologue_ID.split(delim)
                else :
                    data2  = ["NA"]*5                         #less interesting-where a protein is related to itself.
                score  = data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]
                cdict[idx]     = {
                "cluster member ID"         : data1[2],
                "cluster member annotation" : data1[3],
                "cluster member startstop"  : data1[4],
                "homologue organism"        : data2[0],
                "homologue contig"          : data2[1],
                "homologue proteinID"       : data2[2],
                "homologue annot"           : data2[3],
                "homologue startstop"       : data2[4],
                "BLAST bitscore"            : score,
                "CDD domains"               : ",".join(annot_list),
                "Domain definitions"        : ",".join(annot_def_list)  }
                idx += 1
        cframe = pd.DataFrame(cdict).transpose()
        cframe.sort_values('cluster member ID', ascending = False, inplace = True)
        cframe.to_csv(
        abspath(tab_output_dir + "/" + species + "-" + WP_no + ".tsv"),
                sep = "\t")
