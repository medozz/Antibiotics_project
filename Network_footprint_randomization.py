import sys
import igraph
import random as rnd

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
    while i < int(number_of_randomization):
        number_of_seeds = random_list[rnd.randint(0, (len(random_list)-1))]
        number_of_nodes = len(graph.vs)
        seed = 0
        start_node_list = []
        while seed < number_of_seeds:
            node_id = rnd.randint(0, number_of_nodes)
            if node_id not in start_node_list:
                start_node_list.append(node_id)
                seed += 1

        i += 1
