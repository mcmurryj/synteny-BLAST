import networkx as nx
import os
def write_networks(data_box, g_output_dir) :
	cluNet = nx.Graph()
	protNet = nx.Graph()

	for cluster_ID in data_box.keys():
		cluNet.add_node(cluster_ID)

	for cluster_ID in data_box.keys():
		for cluster_member_ID in data_box[cluster_ID].keys():
			for homologue_ID in data_box[cluster_ID][cluster_member_ID].keys():
				score               = float(data_box[cluster_ID][cluster_member_ID][homologue_ID]["bitscore"])
				homologue_parent_ID = data_box[cluster_ID][cluster_member_ID][homologue_ID]["homologue_parent"]
				IDtup               = tuple(sorted([cluster_ID, homologue_parent_ID]))
				IDtup2              = tuple(sorted([cluster_member_ID, homologue_ID]))
				protNet.add_node(IDtup2[0])
				protNet.add_node(IDtup2[1])
				protNet.add_edge(IDtup2[0], IDtup2[1], weight = score)
				if IDtup in cluNet.edges():
					cluNet[IDtup[0]][IDtup[1]]["weight"] += score
				else:
					cluNet.add_edge(IDtup[0], IDtup[1], weight = score)

	#tidy up the network by removing self-connections and deleting loner nodes.
	for n in protNet.nodes():
		if len(protNet[n].keys()) <= 1:
			protNet.remove_node(n)
		if protNet.has_edge(n, n):
			protNet.remove_edge(n,n)

	for n in cluNet.nodes():
		if cluNet.has_edge(n, n):
			cluNet.remove_edge(n,n)

	nx.write_gml(cluNet, os.path.abspath(g_output_dir + "/"+ "graph_o_clusters.gml"))
	nx.write_gml(protNet, os.path.abspath(g_output_dir + "/"+ "graph_o_proteins.gml"))
