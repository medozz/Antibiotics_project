"""
This script will construct the interaction networks using the existing resources. The possible
mapping issues will be saved separately.
Used synonyms will be uniprot IDs
First try to use object orientated programming :)

"""

import re
def collect_targeted_genes_from_Oz_et_al(targeted_genes_link_containing_file):
    """
    This script will give us the BioCyc identifier fo each gene in Oz et al
    :param targeted_genes_link_containing_file: The containing file
    :return: list of ids
    """
    inp = open(targeted_genes_link_containing_file, "rb")
    inp.readline() #header
    eglist = []
    for line in inp:
        line=line.strip()
        egid = re.findall(r"=(EG[0-9]*)", line)
        if egid:
            eglist.append(egid[0])
    return eglist

#elist = collect_targeted_genes_from_Oz_et_al("links.txt")


class edgelist:
    """
    This class contains and edgelist object which is undirected representation of the graph
    """
    def __init__(self):
        self.edges = set()
        self.nodes = set()
        self.properties = {}

    def create_checlist(self, file_):
        self.chcklist_nodes = set()
        inp = open(file_)
        inp.readline()
        for line in inp:
            line = line.split("\t")
            self.chcklist_nodes.add(line[0])
        return self

    def unite_edgelist(self, edgelist):
        """
        Unites two edgelist object
        :return: united edgelist
        """
        self.nodes.union(edgelist.nodes)
        self.edges.union(edgelist.edges)
        return self

    def add_edge(self,add1,add2):
        """
        Makes an edgelist for graph construction. It filrters out self lops and direction of the edges
        :param  self edgelist
        :param add1: node1
        :param add2: node2
        :return: edgelist
        """
        if add1 in self.chcklist_nodes and add2 in self.chcklist_nodes:
            self.nodes.add(add1)
            self.nodes.add(add2)
            if (add1 + " " + add2 not in self.edges) and (add2 + " " + add1 not in self.edges) and (add1 != add2):
                self.edges.add(add1 + " " + add2)
        return self

    def add_property(self, property_name, **kwargs):
        """
        Creates a property for the existing edgelist to fulfill.
        :param property_name: The name of the property
        :param kwargs: variable_type = string(1), int(0) or float(2)
        variable= the existing default variable all the edges will have it.
        :return: self
        """
        self.properties[property_name]={}
        check = 0
        if "variable_type" in kwargs:
            if kwargs["variable_type"] == "string" or kwargs["variable_type"] == 1:
                check = 1
            if kwargs["variable_type"] == "int" or kwargs["variable_type"] == 0:
                check = 0
            if kwargs["variable_type"] == "float" or kwargs["variable_type"] == 2:
                check = 2

        if "variable" in kwargs:
            for edge in self.edges:
                self.properties[property_name][edge] = kwargs["variable"]

        for edge in self.edges:
            if check == 0:
                self.properties[property_name][edge] = 0
            if check == 1:
                self.properties[property_name][edge] = ""
            if check == 2:
                self.properties[property_name][edge] = 0.0
        return self

    def add_property_value(self, property_name,edge_ID,value):
        """
        Add a value to the existing property. The property have to be constructed with the
        add_property function.
        :param property_name: name of the property
        :param edge_ID: the ID of the edge from the edgelist
        :param value: the value of the property
        :return: self
        """
        self.properties[property_name][edge_ID] = value
        return self

    def write_ncol(self, filename, *args):
        out = open(filename, "wb")
        for edge in self.edges:
            out.write(edge)
            if args:
                if type (args) != "string":
                    for property in args:
                        out.write(" "+str(self.properties[property][edge]))
                else:
                    out.write(" " + str(self.properties[property][edge]))
            out.write("\n")
        out.close()








def import_interactions_form_Rajagopola_et_al(file_name, edgelist):
    """
    This function will import the Rajagopola et al file. It will use the text file and construct
    an edge list object. In the edge list object the nodeas are spearated with space and the nodes are
    uniprot IDs. The edge list object has three arguments edges, nodes, properties. Properties contains edge
    properties.There are no duplications in any dirrection, neither self loops.
    :param file_name: the Rajagopala et al file name constructed from Access
    :param edgelist: the original edge list in an edgelist object
    :return: edgelist
    """
    inp = open(file_name)
    inp.readline() #Header
    for line in inp:
        line=line.strip()
        line=line.split("\t")

        if line[0] != "NoData" and line[2] != "NoData":
            edgelist.add_edge(line[0], line[2])
    return edgelist


def search_uniprot_id(inputstring):
    ids = re.findall(r"uniprotkb:(\S{6})", inputstring)
    if ids:
        uid = ids[0]
    else:
        uid = 0
    return uid

def import_interactions_form_DIP(dip_file, edgelist):
    """
    This function will import the DIP E. coli interaction file. It will use the text file and add
    edges to the edglist object. It uses only uniprot IDs. There are no duplications in any
    dirrection, neither self loops.
    :param dip_file: the dip file
    :param edgelist: the edglist object
    :return: edgelist
    """
    inp = open(dip_file)
    inp.readline() #header
    for line in inp:
        line=line.strip()
        line = line.split("\t")
        uniprot_ID1 = search_uniprot_id(line[0])
        uniprot_ID2 = search_uniprot_id(line[1])
        if uniprot_ID1 and uniprot_ID2:
            edgelist.add_edge(uniprot_ID1, uniprot_ID2)
    return edgelist

def import_uniprot_transltaion_dictionarry(uniprot_file):
    inp = open(uniprot_file)
    inp.readline()
    uniprot_dictionarry = {}
    for line in inp:
        line = line.strip()
        line = line.split("\t")
        if line[0] not in uniprot_dictionarry:
            uniprot_dictionarry[line[0]] = [line[2]]
        else:
            uniprot_dictionarry[line[0]].join(line[2])
    return uniprot_dictionarry

def import_string_edges(stringfile, edgelist, uniprot_translation_dictionarry):
    inp = open(stringfile)
    for line in inp:
        line = line.strip()
        line = line.split(" ")
        if float(line[6])>0 or float(line[7])>0:
            if line[0].upper() in uniprot_translation_dictionarry and line[1].upper() in uniprot_translation_dictionarry:
                for uniprot_ID1 in uniprot_translation_dictionarry[line[0].upper()]:
                    for uniprot_ID2 in uniprot_translation_dictionarry[line[1].upper()]:
                        #print uniprot_ID2, uniprot_ID1
                        edgelist.add_edge(uniprot_ID1, uniprot_ID2)
    return edgelist




uniprot_translation_dictionarry =import_uniprot_transltaion_dictionarry("STRING_to_uniprot.tab")

elist = edgelist() #constructor
elist.create_checlist("uniprot-e_coli_k12_07_11_2016.tab")
elist = import_string_edges("511145string_data.txt", elist ,uniprot_translation_dictionarry)
print "Number of edges after STRING:", len(elist.edges)
print "Number of nodes after STRING:", len(elist.nodes)

elist = import_interactions_form_Rajagopola_et_al("usable_interactions.txt", elist)
print "Number of edges after Rajagopola et al:", len(elist.edges)
print "Number of nodes after Rajagopola et al:", len(elist.nodes)

elist = import_interactions_form_DIP("Ecoli20160731.txt", elist)
print "Number of edges after DIP:", len(elist.edges)
print "Number of nodes after DIP", len(elist.nodes)


elist = import_interactions_form_DIP("E_coli_MPIDB_export_12_10_16.tab", elist)
print "Number of edges after MPIDB:", len(elist.edges)
print "Number of nodes after MPIDB:", len(elist.nodes)
elist.write_ncol("30_01_2017_DIP_RAJAGOPOLA_MPIDB_STRING_E_coli_K12_PPI.ncol")
