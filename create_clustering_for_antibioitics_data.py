from pandas import DataFrame, Series
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from sklearn import cluster, datasets
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster.bicluster import SpectralBiclustering
import mca
import igraph


from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import distance

def coloring_dendogram(dataset, list_of_ids, type_of_data, instances_colors_dictionarry):
    """
    :param dataset: imported dataset from a csv file
    :param list_of_ids: id list (patient IDs), leaves of the clustering diagram
    :param type_of_data: The type of data "continues" if values, "discrete" if discrete instances
    :param instances_colors_dictionarry: If type_of_data continues two RGB tuple for min and max. If type_of_data
            "discrete" then a dictionarry for each instances in data and an RGB tupple
    :return: a list of IDS and colors dictionarry
    """
    label_colors = {}
    values = []
    for patient in list_of_ids:
        print patient
        value_row =  [dataset["ID"] == patient]
        print value_row
        value = value_row["Value"].tolist()
        if type_of_data == "discrete":
            for instance in instances_colors_dictionarry:
                print value
                if value[0] == instance:
                    label_colors[patient] = instances_colors_dictionarry[instance]

        if type_of_data == "continues":
            values.append(float(value[0]))

    if type_of_data == "continues":
        minimum = min(values)
        maximum = max(values)
        print minimum
        dr = float(instances_colors_dictionarry[1][0]-instances_colors_dictionarry[0][0])
        dg = float(instances_colors_dictionarry[1][1]-instances_colors_dictionarry[0][1])
        db = float(instances_colors_dictionarry[1][2]-instances_colors_dictionarry[0][2])

        rmin = float(instances_colors_dictionarry[0][0])
        gmin = float(instances_colors_dictionarry[0][1])
        bmin = float(instances_colors_dictionarry[0][2])

        dv = maximum-minimum
        k = 0
        while k <len(values):
            value = values[k]
            r = rmin + ((value - minimum) / dv) * dr
            g = gmin + ((value - minimum) / dv) * dg
            b = bmin + ((value - minimum) / dv) * db
            label_colors[list_of_ids[k]]=(r,g,b)
            k = k+1

    return label_colors


def make_dendogram(panda_df, file_name, label_colors):
    """
    It makes a dendogram using Hamming distance from the previously setted data using avarage method.
    :param panda_df: the source data frame
    :param file_name: the output file name
    :param label_colors: the coloring option of the label in current form red and blue
    """
    #print panda_df.columns.tolist()
    data = panda_df.values
    data = data.T
    Z = linkage(data, method="average",  metric="hamming") #metric='jaccard') #euclidean')#"jaccard"
    Y = distance.pdist(data, metric="hamming") #metric='jaccard') # metric="hamming") #'euclidean'
    W = distance.pdist(data.T, metric="jaccard")
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical clustering of antibiotics according to affected nodes in the E. coli K12 PPI network')
    plt.xlabel('Antibioitics')
    #plt.ylabel('Jaccard Distance')
    #plt.ylabel('Euclidean Distance')
    plt.ylabel('Hamming Distance')
    #print gender_dataset["PatientID"]

    print label_colors
    if not label_colors:
        label_colors = {}
        labels = panda_df.columns.tolist()
        for label in labels:
            print label
            label_colors[label] = (0,0,0)

    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.0,  # font size for the x axis labels
        color_threshold= 0.18, #6, 0.10, #4.5, #150, #0.6,
        labels=panda_df.columns.tolist()
    )
    print Z[0]
    #plt.show()
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(label_colors[lbl.get_text()])
    plt.savefig(file_name, dpi =600, format="png", transparent=True)

    k = 6

    clusters = fcluster(Z, k, criterion='maxclust')
    out = open(file_name.replace(".png","_cluster.txt"),"wb")
    out2 = open(file_name.replace(".png", "_distance_matrix.txt"), "wb")
    Y = distance.squareform(Y)
    for k in range(len(clusters)):
        out.write(str(clusters[k])+"\t"+str(panda_df.columns.tolist()[k])+"\n")
        for i in range(k):
            out2.write(str(panda_df.columns.tolist()[k])+"\t"+str(panda_df.columns.tolist()[i])+"\t")
            out2.write(str(1-Y[i,k]))
            out2.write("\n")

    out3 = open(file_name.replace(".png", "protein_distance_matrix.txt"), "wb")
    W = distance.squareform(W)
    for k in range(len(panda_df.index.tolist())):
        for i in range(k):
            out3.write(str(panda_df.index.tolist()[k]) + "\t" + str(panda_df.index.tolist()[i]) + "\t")
            out3.write(str(1.0-(W[i, k])))
            out3.write("\n")
    out.close()
    out2.close()
    out3.close()
def create_first_neighbours (G, antibioitics_dataset, out):

    antibioitics_set = set()
    proteins_set = set()
    for antibiotics in antibioitics_dataset.columns:
        G.vs[antibiotics] = 0
        antibioitics_set.add(antibiotics)
        for protein in antibioitics_dataset.index:
            if antibioitics_dataset.at[protein, antibiotics] > 0:
                if protein in G.vs["name"]:
                    node = G.vs.select(name_eq=protein)
                    node[antibiotics] = 1
                    print node["name"]
                    proteins_set.add(node["name"][0])
                    neighbours = G.neighborhood(node[0])
                    # This part requires a lot to do. It is unknown what information flow method would be the best.
                    # The current method is bit based, if a protein is affected then it is affected whatever way.
                    for neighbour in neighbours:
                        print neighbour
                        if G.vs[neighbour][antibiotics] != 1: #I know it is useless for future.
                            G.vs[neighbour][antibiotics] = 1
                            proteins_set.add(G.vs[neighbour]["name"])

    out = open("antibioitics_resultfile.txt", "wb")
    proteins_list = list(proteins_set)
    antibioitics_list = list(antibioitics_set)
    out.write("\t"+"\t".join(antibioitics_list)+"\n")
    for protein in proteins_list:
        out.write(protein)
        node = G.vs.select(name_eq=protein)
        for antibiotics in antibioitics_list:
            out.write("\t"+str(node[antibiotics][0]))
        out.write("\n")
    out.close()

G = igraph.Graph.Read_Ncol("07_11_2016_DIP_RAJAGOPOLA_MPIDB_E_coli_K12_PPI.ncol")
antibioitics_dataset = pd.read_csv("Antibioitics_names_distribution.txt", index_col=0, sep="\t")
out = "antibioitics_resultfile.txt"
#create_first_neighbours(G, antibioitics_dataset, out)

antibioitics_dataset = pd.read_csv("antibioitics_resultfile.txt", index_col=0, sep="\t", header=0)
make_dendogram(antibioitics_dataset, "antibiotics_clusters_Hamming_2.png", 0)







