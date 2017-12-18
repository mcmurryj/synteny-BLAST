#!/usr/bin/env python
'''Utilities for plotting gene clusters.
   Todo: integrate into main pipeline.
   Todo: Fix genbank skinner to preserve sense info.
   Todo: Flop things around so that the cluster_ID gene is always on the + strand;
   this will make it easier to visualize synteny.
   Todo: Fix stuff so that identical cluster_member_IDs are rejected.'''
import operator
from os.path import abspath
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

#check http://biopython-cn.readthedocs.io/zh_CN/latest/en/chr17.html
#http://biopython.org/DIST/docs/GenomeDiagram/userguide.pdf
#create a diagram object,  using its methods add track(s), and use the track methods to add feature-sets

def getClustSS(clusterMembers):
    cstart = 9E10
    cstop = 0
    for cm in clusterMembers:
        startstop = cm.split("|")[4]
        startstop = re.sub('[<>$%^&*#|\[\]]', "", startstop)
        start = int(startstop.split("-")[0])
        stop = int(startstop.split("-")[1])
        if(start < cstart):
            cstart = start
        if(stop > cstop):
            cstop = stop
    output = (cstart, cstop, cstop-cstart)
    return output

def getMaxSpan(sortedClusters):
    """Function to get the maximum span from a list of
       cluster_IDs, in terms of base pair length of cluster."""
    maxSpan = 0
    for e in sortedClusters:
        clusterMembers = data_box[e].keys()
        span = getClustSS(clusterMembers)[2]
        if span > maxSpan:
            maxSpan = span
    return maxSpan

def getTopClusters(cluNet, node, ncutoff=8):
    """Function that returns the clusters that are most closely related
    to the cluster of interest.
    Input:
    NetworkxClusterGraph[node]
    Output:
    List of cluster IDs.
    TODO:  remove clusters keyed on identical protein IDs.
    """
    neigh = cluNet[node]
    #make dict of weights keyed on neighbor IDs.
    edgeWeights = {n:neigh[n]["weight"] for n in neigh}
    #Add the node its own self with superhigh value weights
    edgeWeights[node] = 999999
    #Sort edgeWeights.items on the second element of each item, highest values descending to lowest
    sortedEdges = sorted(edgeWeights.items(), key=operator.itemgetter(1), reverse=True)
    #go back throgh and remove clusters keyed on identical proteins.
    resortedEdges = []
    lastCMID = ""
    for se in sortedEdges:
        wpID = se[0].split("|")[2]
        if wpID != lastCMID:
            resortedEdges.append
    sortedClusters = tuple([e[0] for e in sortedEdges])[0:ncutoff]


    return (sortedClusters, sortedEdges)

def printClusterPics(cluNet, data_box, annot_def_dict, output_dir, ncutoff=8):
    """This abomination takes in a network of gene clusters, a dict of cluster_member_ID:pfam domains,
       and a place to dump output.  It goes thru and prints pictures of each cluster and its best buddies.
       Possible defects:
       Proteins can be related, but won't be colored if there is no pfam domain annotation.
       Handling of multidomain proteins is very crude.
       Also, much much too long."""
    for node in cluNet.nodes():
        #record the species for file printing poirposes
        parentSpecies, parentAcession, parentWpID, parentAnnotation, parentStartstop = node.split("|")
        #Make a new diagram.  TODO add the path for the output and whatnot.
        gd_diagram = GenomeDiagram.Diagram()
        #Return the top most related clusters.
        sortedClusters = getTopClusters(cluNet, node)[0]
        #the largest nucleotide span in all the clusters under consideration.
        maxSpan = getMaxSpan(sortedClusters)
        #Get info about pfam domains in all the clusters under consideration
        id_list = makeIDList(sortedClusters, data_box)
        domainCounts = getPfamCounts(annot_def_dict, id_list)
        colorDict = getColorDict(domainCounts)
        # for each of the best clusters:
        for e in sortedClusters:
            #Make the track and add a feature set to the track. scale = 0 turns off the scale, greytrack adds a background.  What is the 1 for?  Who knows.
            gd_track_for_features = gd_diagram.new_track(1, name="", scale=0, greytrack=True)
            gd_feature_set = gd_track_for_features.new_set()
            #get the cluster_member_IDs from data_box using the cluster_ID
            #TODO this will be replaced with a clusterID:clusterMemberID dict based on the BLAST2 db? maybe.
            clusterMembers = data_box[e].keys()
            #start stop span
            sss = getClustSS(clusterMembers)
            offset = int((maxSpan - sss[2])/2 -sss[0])
            for cm in clusterMembers:
                organism, acession, wpID, annotation, startstop = cm.split("|")
                startstop = re.sub('[<>$%^&*#|\[\]]', "", startstop)
                q_start = int(startstop.split("-")[0])
                q_stop = int(startstop.split("-")[1])
                #encode the directional info.
                if q_start < q_stop:
                    q_strand = 1
                    label_angle = 15
                elif q_start > q_stop:
                    q_strand = -1
                    label_angle = 179.9
                #else:
                #    strand = None
                #Uh, I think this should make the whole shebang roughly centered?
                feature = SeqFeature(FeatureLocation(q_start + offset, q_stop + offset),
                                     ref=wpID, strand=q_strand)
                #What color should it be?
                pfam = domainHash(annot_def_dict[cm])
                if len(pfam) > 0 and pfam in colorDict:
                    color = colorDict[pfam]    #pfam is a tuple of domain names, usually just one
                else:
                    color = "0x000000"  #paint it black
                gd_feature_set.add_feature(feature, sigil="ARROW", color=color,
                                           label=True, label_position="start",
                                           name=feature.ref, label_size = 7,
                                           label_angle=label_angle)
        #TODO Add a legend at the bottom, too. each pfam dom is a colored square. ref = name of pfam domain.
        gd_track_for_features = gd_diagram.new_track(0, name="", scale=0, greytrack=False)
        gd_feature_set = gd_track_for_features.new_set()
        counter = 0
        increment = maxSpan/len(colorDict.keys())
        for domain_tup in colorDict:
            q_start = round(counter * increment)
            q_stop =  round((counter + .33) * increment)
            counter += 1
            feature = SeqFeature(FeatureLocation(q_start, q_stop),
                                 ref=domain_tup[0], strand=None)
            gd_feature_set.add_feature(feature, sigil="BOX", color=colorDict[domain_tup],
                                       label=True, name=feature.ref,
                                       label_size = 7, label_angle=15)
        #Done reading in all the info! Now write the plot
        name = abspath(output_dir + "/" + parentSpecies + "-" + parentWpID + ".pdf")
        drawMap(gd_diagram, name, plotWidth=maxSpan)


def drawMap(gd_diagram, plotName, plotWidth=0):
    gd_diagram.draw(format="linear",
                    pagesize='A4',
                    fragments=1,
                    start=0, end=plotWidth,
                    track_size=.5)
    gd_diagram.write(plotName, "PDF")


def makeIDList(sortedClusters, data_box):
    """Take a list of Cluster IDs.
       Return a list of all cluster member IDS in those clusters"""
    alist = []
    for c in sortedClusters:
        alist += list(data_box[c].keys())
    aset = set(alist)
    return list(aset)


def getPfamCounts(annot_def_dict, id_list):
    """Accepts dictionary of protein ids:domain annotations,
       and list of protein ids.  Returns a tuple of domain annotation-count pairs,
       ordered on the number of counts.
       May want to reconsider treatment of multi-domain proteins. """
    pfamCounts = {}
    #Iterate thru the list of protein IDs
    for id in id_list:
        #If the id is annotated
        if id in annot_def_dict and len(annot_def_dict[id]) > 0:
            pfam = domainHash(annot_def_dict[id])
            if pfam == ("none") :
                continue
            if pfam not in pfamCounts:
                pfamCounts[pfam] = 1
            elif pfam in pfamCounts:
                pfamCounts[pfam] += 1
    sortedPfams = sorted(pfamCounts.items(), key=operator.itemgetter(1), reverse=True)
    return sortedPfams

def domainHash(pfam_list):
    if len(pfam_list) == 0:
        return ("none")
    else:
        return tuple(sorted(pfam_list))

def getColorDict(domainCounts):
    '''Horribly unstylish function to map pfam domains to colors.
    Hex code color sets that are relatively distinguishable by most colorblind people.
    Color schemes from Paul Tol: https://personal.sron.nl/~pault/'''
    ptolthree = ['0x4477AA', '0xDDCC77', '0xCC6677']
    ptolfour = ['0x4477AA', '0x117733', '0xDDCC77', '0xCC6677']
    ptolfive = ['0x332288', '0x88CCEE', '0x117733', '0xDDCC77', '0xCC6677']
    ptolsix = ['0x332288', '0x88CCEE', '0x117733', '0xDDCC77', '0xCC6677', '0xAA4499']
    ptolseven  = ['0x332288', '0x88CCEE', '0x44AA99',
                    '0x117733', '0xDDCC77', '0xCC6677', '0xAA4499']
    ptoleight = ['0x332288', '0x88CCEE', '0x44AA99', '0x117733',
                    '0x999933', '0xDDCC77', '0xCC6677', '0xAA4499']
    ptolnine = ['0x332288', '0x88CCEE', '0x44AA99', '0x117733',
                    '0x999933', '0xDDCC77', '0xCC6677', '0x882255', '0xAA4499']
    ptolten = ['0x332288', '0x88CCEE', '0x44AA99', '0x117733', '0x999933',
                  '0xDDCC77', '0x661100', '0xCC6677', '0x882255', '0xAA4499']
    ptoleleven = ['0x332288', '0x6699CC', '0x88CCEE', '0x44AA99', '0x117733', '0x999933',
                     '0xDDCC77', '0x661100', '0xCC6677', '0x882255', '0xAA4499']
    if len(domainCounts) >= 11:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptoleleven)}
    if len(domainCounts) == 10:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolten)}
    if len(domainCounts) == 9:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolnine)}
    if len(domainCounts) == 8:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptoleight)}
    if len(domainCounts) == 7:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolseven)}
    if len(domainCounts) == 6:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolsix)}
    if len(domainCounts) == 5:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolfive)}
    if len(domainCounts) == 4:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolfour)}
    if len(domainCounts) == 3:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolthree)}
    if len(domainCounts) == 2:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolthree[0:2])}
    if len(domainCounts) == 1:
        col_dict = {id[0]:col for id, col in
                      zip(domainCounts, ptolthree[0])}

    return col_dict
