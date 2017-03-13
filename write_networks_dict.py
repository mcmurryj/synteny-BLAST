import networkx as nx
import os
def write_networks(data_box, g_output_dir, annot_dict = False) :
    cluNet  = {}
    protNet = {}
    for cluster_ID in data_box.keys():
        cluNet[cluster_ID] = {}


    for cluster_ID in data_box.keys():
        for cluster_member_ID in sorted(data_box[cluster_ID].keys()):
            for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
                #Important-score must be a float, otherwise adding gets funky
                score               = float(data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"])
                homologue_parent_ID = data_box[cluster_ID][cluster_member_ID][homologue_ID]["homologue_parent"]
                IDtup               = tuple(sorted([cluster_ID, homologue_parent_ID]))
                IDtup2              = tuple(sorted([cluster_member_ID, homologue_ID]))
                protNet[IDtup2[0]]  = {}
                protNet[IDtup2[1]]  = {}
                protNet[IDtup2[0]]  = {IDtup2[1] : {"weight" :  score}}
                if IDtup[1] in cluNet[IDtup[0]].keys():
                    cluNet[IDtup[0]][IDtup[1]]["weight"]     += score
                else:
                    cluNet[IDtup[0]][IDtup[1]] =   {"weight" :  score}

    cluG = nx.Graph(cluNet)
    for n in cluG.nodes():
        if cluG.has_edge(n, n):
            cluG.remove_edge(n,n)

    protG = nx.Graph(protNet)
    for n in protG.nodes() :
        if protG.has_edge(n, n) :
            protG.remove_edge(n, n)

    nx.write_gml(cluG,  os.path.abspath(g_output_dir + "/"+ "graph_o_clusters.gml"))
    nx.write_gml(protG, os.path.abspath(g_output_dir + "/"+ "graph_o_proteins.gml"))
