"""Functions for calculating and printing summary data."""
from itertools import product
from os.path import abspath
from glob import glob
import pandas as pd

def summarize(tab_output_dir) :
    """Read the raw output, then calculate and print summary data."""
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
            tabsum.sort_values('No._of_homologues ', ascending = False, inplace = True)
            tabsum.index.rename("NCBI protein ID", inplace = True)
            tabsum.to_csv(abspath(f[:-4] + "synop.tsv"), sep = "\t")

#def print_tables(data_box, tab_output_dir, annot_dict = {}, annot_def_dict = {}):
#    """Iterate thru all the info in data_box and print to tab_output_dir"""
#    delim  = "|"
#    for cluster_ID in data_box.keys():
#        cdict         = {}
#        #get data chunks from the titles of the fasta sequences
#        data0         = cluster_ID.split(delim)
#        species       = data0[0].replace(" ", "_")
#        WP_no         = data0[2]                #third element is the WP_number
#        org_cont      = data0[0:2]                #first and second elements are the species, contig
#        idx           = 0
#        #Loop thru the dict and print out the good stuff.
#        for cluster_member_ID in data_box[cluster_ID].keys():
#            data1 = cluster_member_ID.split(delim)
#            #Check to see if there is a domain annotation list; if there is add it.
#            if cluster_member_ID in annot_dict.keys() :
#                annot_list     = annot_dict[cluster_member_ID]
#                annot_def_list = annot_def_dict[cluster_member_ID]
#            else :
#                annot_list     = []
#                annot_def_list = []
#            for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
#                if cluster_member_ID != homologue_ID :       #if it is not self-identity
#                    data2 = homologue_ID.split(delim)
#                else :
#                    data2  = ["NA"]*5                         #less interesting-where a protein is related to itself.
#                score  = data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"]
#                cdict[idx]     = {
#                "cluster member ID"         : data1[2],
#                "cluster member annotation" : data1[3],
#                "cluster member startstop"  : data1[4],
#                "homologue organism"        : data2[0],
#                "homologue contig"          : data2[1],
#                "homologue proteinID"       : data2[2],
#                "homologue annot"           : data2[3],
#                "homologue startstop"       : data2[4],
#                "BLAST bitscore"            : score,
#                "CDD domains"               : ",".join(annot_list),
#                "Domain definitions"        : ",".join(annot_def_list)  }
#                idx += 1
#        #Convert the dictionary with data at the homologue level into a dataframe
#        cframe = pd.DataFrame(cdict).transpose()
#        #Sort the data at the level of cluster member
#        cframe.sort_values('cluster member ID', ascending = False, inplace = True)
#        #Write the data to csv file
#        cframe.to_csv(
#        abspath(tab_output_dir + "/" +
#                species + "-" +
#                WP_no + ".tsv"),
#                sep = "\t")

def contingencyTable(tab_output_dir):
    """Idea is to print contingency tables for PFAM domain co-incidence.
       Use synopsis tab output."""
    #Make a blank dictionary
    contingency = {}
    #synopsis files
    synop_files = glob(abspath(tab_output_dir + "/*synop.tsv"))
    for sf in synop_files:
        #Read the file
        syn_tab = pd.read_csv(sf, sep="\t", header=0)
        #Blank the domain list- we are switching to a different cluster.
        domain_list = []
        for domains in syn_tab['Pfam_domain_def.']:
            domain_list += str(domains).split(",")
        #Summarize cluster-level results.
        #For every pair of PFAM domains in the list of annotations....
        for pA, pB in product(domain_list, repeat=2):
            #If we've never seen domain A before
            if not pA in contingency.keys():
                #Make a blank entry
                contingency[pA] = {}
            #if we've never seen the pA->pB combo...
            if not pB in contingency[pA].keys():
                #count is now 1
                contingency[pA][pB] = 1
            else:
                #If we have seen pA->pB before, increment counter 1
                contingency[pA][pB] += 1
    #Convert to table:
    contingency_frame = pd.DataFrame(contingency)
    #Fill NAs with 0
    contingency_frame.fillna(0, inplace=True)
    #margins
    cont_rowsum = contingency_frame.sum(axis=0)
    cont_rowsum.sort_values(ascending=False, inplace=True)
    contingency_frame = contingency_frame.loc[cont_rowsum.index]
    contingency_frame = contingency_frame[cont_rowsum.index]
    #write contingency_frame
    contingency_frame.to_csv(abspath(tab_output_dir + "/" +
                                     "contingency_table.tsv"),
                            sep = "\t")
