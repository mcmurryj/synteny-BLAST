import networkx as nx
import os
def write_networks(data_box, g_output_dir='') :
    cluNet = nx.Graph()
    protNet = nx.Graph()

    for cluster_ID in data_box.keys():
        cluNet.add_node(cluster_ID)

    for cluster_ID in data_box.keys():
        sorted_cmemb = sorted(data_box[cluster_ID].keys())
        for cluster_member_ID in sorted_cmemb:
            for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
                #Important-score must be a float, otherwise adding gets funky
                score               = float(data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"])
                homologue_parent_ID = data_box[cluster_ID][cluster_member_ID][homologue_ID]["homologue_parent"]
                IDtup               = tuple(sorted([cluster_ID, homologue_parent_ID]))
                IDtup2              = tuple(sorted([cluster_member_ID, homologue_ID]))
                protNet.add_node(IDtup2[0])
                protNet.add_node(IDtup2[1])
                protNet.add_edge(IDtup2[0], IDtup2[1], weight = score)
                if IDtup in cluNet.edges():                     #networkx seems to run really really slow...
                    cluNet[IDtup[0]][IDtup[1]]["weight"] += score
                else:
                    cluNet.add_edge(IDtup[0], IDtup[1], weight = score)

    #tidy up the network by removing self-connections and deleting loner nodes.
    for n in protNet.nodes():
        if len(protNet[n].keys()) <= 1:
            protNet.remove_node(n)
        if protNet.has_edge(n, n):
            protNet.remove_edge(n,n)

    #remove self edges.
    for n in cluNet.nodes():
        if cluNet.has_edge(n, n):
            cluNet.remove_edge(n,n)
    #remove weightless edges.
    for e1, e2 in cluNet.edges():
        if cluNet[e1][e2]["weight"] < 1:
            cluNet.remove_edge(e1,e2)

    #If output folder is specified print print the stuff, in any case return the network.
    return cluNet
    if g_output_dir != '':
        nx.write_gml(cluNet,  os.path.abspath(g_output_dir + "/"+ "graph_o_clusters.gml"))
        nx.write_gml(protNet, os.path.abspath(g_output_dir + "/"+ "graph_o_proteins.gml"))


def getNetwork(data_box):
    '''Make networkX network from data_box and return the network.
       Should not be needed, kept only for debugging purposes.'''
    cluNet = nx.Graph()

    for cluster_ID in data_box.keys():
        cluNet.add_node(cluster_ID)

    for cluster_ID in data_box.keys():
        sorted_cmemb = sorted(data_box[cluster_ID].keys())
        for cluster_member_ID in sorted_cmemb:
            for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
                #Important-score must be a float, otherwise adding gets funky
                score               = float(data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"])
                homologue_parent_ID = data_box[cluster_ID][cluster_member_ID][homologue_ID]["homologue_parent"]
                IDtup               = tuple(sorted([cluster_ID, homologue_parent_ID]))
                IDtup2              = tuple(sorted([cluster_member_ID, homologue_ID]))
                if IDtup in cluNet.edges():                     #networkx seems to run really really slow...
                    cluNet[IDtup[0]][IDtup[1]]["weight"] += score
                else:
                    cluNet.add_edge(IDtup[0], IDtup[1], weight = score)

    #remove self edges.
    for n in cluNet.nodes():
        if cluNet.has_edge(n, n):
            cluNet.remove_edge(n,n)
    #remove weightless edges.
    for e1, e2 in cluNet.edges():
        if cluNet[e1][e2]["weight"] < 1:
            cluNet.remove_edge(e1,e2)

    return cluNet
