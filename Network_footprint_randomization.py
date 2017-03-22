# Copyright (c) 2017 Dezso Modos
# University of Cambridge
# Licensed under the terms of the MIT License

import sys
import igraph
import random as rnd
import copy
import pandas as pd
from pandas import DataFrame, Series


if __name__ == '__main__':
    """
    This script is a simple random network footprint creator. It requires an graph as an .ncol file (this is the
    simplest which iGraph can process). This is the first argument. Then it requires a file with number of original
    nodes separated by new lines. Ideally it is a csv or txt extension, but it does not matter.
    The third argument is the randomization number. The output file will concatenate the fieldnames adn will
    write random at the end.
    """
    rnd.seed = "I am the random seed"  # Yes it can be string :)

    graph = sys.argv[1]
    affected_targets_number_distribution = sys.argv[2]
    number_of_randomization = sys.argv[3]
    filename_help = sys.argv[2].replace(".txt","")
    filename_help = filename_help.replace(".csv","")
    filename_help = filename_help.replace(".dat", "")
    outfile_name = sys.argv[1].replace(".ncol","")+"_"+filename_help+"_random.txt"

    try:
        graph = igraph.Graph.Read_Ncol(graph)
        inp = open(affected_targets_number_distribution)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

    random_list = []
    for line in inp:
        try:
            line = int(line.strip())
        except ValueError:
            print "Could not convert all values form", str(sys.argv[2]), " to integers."
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise
        random_list.append(line)

    try:
        number_of_randomization = int(number_of_randomization)
    except ValueError:
        print "The randomization is not an integer."
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

    i = 0
    occurence_dictionarry = {}

    while i < int(number_of_randomization):
        # Creating the dictionarry to save the data. It can be transformed to pandas data frame.
        occurence_dictionarry[i] = {}
        for vertex_name in graph.vs["name"]:
            occurence_dictionarry[i][vertex_name] = 0
        number_of_seeds = random_list[rnd.randint(0, (len(random_list)-1))]
        number_of_nodes = len(graph.vs)-1
        seed = 0
        start_node_list = []


        while seed < number_of_seeds:
            node_id = rnd.randint(1, number_of_nodes)
            if node_id not in start_node_list:
                start_node_list.append(node_id)
                seed += 1
        graph.vs["value"] = 0
        neighbours_set = set()
        for cr in start_node_list:
            graph.vs[cr]["value"] += 2
            neighbours = graph.vs[cr].neighbors(mode="ALL")
            for node in neighbours and node.index not in neighbours_set:
                node["value"] += 1 # this currently will increase the importance of neighbours only once
                neighbours_set.add(node.index)
        for node in graph.vs:
            occurence_dictionarry[i][node["name"]] = node["value"]
        i += 1
    df = DataFrame.from_dict(occurence_dictionarry)
    df.to_csv(outfile_name, sep="\t")