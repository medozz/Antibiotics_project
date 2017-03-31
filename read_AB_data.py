"""
This script will go through the whole antibiotics network pipeline.
1. read in the E. coli network.
2. Calculate each antibiotics footprint in the PPI sense
3. Calculate each antibiotics TF footprint

Results:
Table for each proteins with weights of AB resistance involvement columns antibioitics
Table of each proteins with type list of AB resistance involvement
Table yes or no involvement for each proteins
Similarity between antibiotics using cosine similarity
Similarity between antibiotics using Tanimoto Similarity
"""
import pandas as pd
import igraph
"""
#This part is the network construction
import Constructing_interaction_network as constin
from Constructing_interaction_network import edgelist
elist = edgelist() #constructor
elist.create_checlist("uniprot-e_coli_k12_07_11_2016.tab")
elist = constin.import_interactions_form_Rajagopola_et_al("usable_interactions.txt", elist)
print "Number of edges after Rajagopola et al:", len(elist.edges)
print "Number of nodes after Rajagopola et al:", len(elist.nodes)

elist = constin.import_interactions_form_DIP("Ecoli20160731.txt", elist)
print "Number of edges after DIP:", len(elist.edges)
print "Number of nodes after DIP", len(elist.nodes)


elist = constin.import_interactions_form_DIP("E_coli_MPIDB_export_12_10_16.tab", elist)
print "Number of edges after MPIDB:", len(elist.edges)
print "Number of nodes after MPIDB:", len(elist.nodes)
elist.write_ncol("30_01_2017_DIP_RAJAGOPOLA_MPIDB_STRING_E_coli_K12_PPI.ncol")
"""

def create_first_neighbours (G, antibioitics_dataset):
    """
    This part will create the first neighbours of the PPi datest adn construct the needed addtions

    :param G: igraph graph
    :param antibioitics_dataset: the involved antibioitics
    :return: the G graph with the prper modifications and a protein_set which contains all affected proteins
    """

    antibioitics_set = set()
    proteins_set = set()
    for antibiotics in antibioitics_dataset.columns:
        for node in G.vs:
            node[antibiotics] = {"types" : set(), "sumvalue" : 0.0, "affected" : 0}
        #print G.vs[antibiotics]
        antibioitics_set.add(antibiotics)
        AB_proteins = set()
        Neighbour_proteins = set()
        for protein in antibioitics_dataset.index:
            if antibioitics_dataset.at[protein, antibiotics] > 0:
                if protein in G.vs["name"] and protein not in AB_proteins:
                    node = G.vs.select(name_eq=protein)
                    node = node[0]
                    #print node
                    node[antibiotics]["types"].add("Mut")
                    node[antibiotics]["sumvalue"] += 2
                    node[antibiotics]["affected"] = 1
                    proteins_set.add(node["name"])
                    AB_proteins.add(node["name"])
                    neighbours = G.neighborhood(node)
                    # This part requires a lot to do. It is unknown what information flow method would be the best.
                    # The current method is bit based, if a protein is affected then it is affected whatever way.
                    for neighbour in neighbours:
                        if G.vs[neighbour]["name"] not in Neighbour_proteins:
                            print G.vs[neighbour]["name"]
                            G.vs[neighbour][antibiotics]["affected"] = 1
                            G.vs[neighbour][antibiotics]["sumvalue"] += 1
                            G.vs[neighbour][antibiotics]["types"].add("FN")
                            proteins_set.add(G.vs[neighbour]["name"])
                            Neighbour_proteins.add(G.vs[neighbour]["name"])

    antibioitics_list = list(antibioitics_set)
    return G,proteins_set,antibioitics_list, Neighbour_proteins


def read_in_tf_data_from_file_to_igraph_graph(file_, G, proteins_set, antibioitics_list):
    """
    This function reads in the transcription factors data as well and add to the graph.
    :param file_: the trasncription data file
    :param G: the igraph graph the PPi data are read in.
    :return: The igrapgh graph.
    """
    inp = open(file_)
    antibiotics_dictionarry_for_TFs = {}
    inp.readline() #header
    for line in inp:
        line = line.strip()
        line = line.split("\t")
        if line[2] not in antibiotics_dictionarry_for_TFs and line[2] in antibioitics_list:
            antibiotics_dictionarry_for_TFs[line[2]] ={}
            antibiotics_dictionarry_for_TFs[line[2]]["TF"] = set()
            antibiotics_dictionarry_for_TFs[line[2]]["target"] = set()
        antibiotics_dictionarry_for_TFs[line[2]]["TF"].add(line[0])
        antibiotics_dictionarry_for_TFs[line[2]]["target"].add(line[1])

    for antibiotic in antibiotics_dictionarry_for_TFs:
        for TF in antibiotics_dictionarry_for_TFs[antibiotic]["TF"]:
            if TF in G.vs["name"]:
                node = G.vs.select(name_eq=TF)
                node = node[0]
                node[antibiotic]["types"].add("TF")
                node[antibiotic]["sumvalue"] += 1
                node[antibiotic]["affected"] = 1
                proteins_set.add(node["name"])
        for target in antibiotics_dictionarry_for_TFs[antibiotic]["target"]:
            if target in G.vs["name"]:
                node = G.vs.select(name_eq=target)
                node = node[0]
                node[antibiotic]["types"].add("target")
                node[antibiotic]["sumvalue"] += 1
                node[antibiotic]["affected"] = 1
                proteins_set.add(node["name"])
                neighbours = G.neighborhood(node)
                for neighbour in neighbours:
                    if "FN" not in G.vs[neighbour][antibiotic]["types"]:
                        G.vs[neighbour][antibiotic]["affected"] = 1
                        G.vs[neighbour][antibiotic]["sumvalue"] += 1
                        G.vs[neighbour][antibiotic]["types"].add("FN")
                        proteins_set.add(G.vs[neighbour]["name"])
    return G, proteins_set


def write_out_graph_dataset(G, proteins_set, antibioitics_list, out_file):

    out = open(out_file, "wb")
    proteins_list = list(proteins_set)
    out.write("\t"+"\t".join(antibioitics_list)+"\n")
    for protein in proteins_list:
        out.write(protein)
        node = G.vs.select(name_eq=protein)
        node = node[0]
        for antibiotics in antibioitics_list:
            out.write("\t"+str(",".join(node[antibiotics]["types"])))
        out.write("\n")
    out.close()

antibioitics_dataset = pd.read_csv("31_03_2017_AB_mutation.txt", index_col=0, sep="\t", header=0)
G = igraph.Graph.Read_Ncol("07_11_2016_DIP_RAJAGOPOLA_MPIDB_E_coli_K12_PPI.ncol")
G,proteins_set,antibioitics_list, Neighbour_proteins = create_first_neighbours(G, antibioitics_dataset)
G, proteins_set = read_in_tf_data_from_file_to_igraph_graph("31_03__2017_TF_targets_datasset.txt",G,proteins_set,antibioitics_list)
write_out_graph_dataset (G,proteins_set,antibioitics_list, "31_03_2017_affected_types.txt")