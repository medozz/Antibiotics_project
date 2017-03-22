import json
import re

CARD = json.load(open("card.json"))
k=0
i=0
Oneset = set()
for id in CARD:
    i=i+1
    try:
        for seqid in CARD[id]["model_sequences"]["sequence"]:
            if CARD[id]["model_sequences"]["sequence"][seqid]['NCBI_taxonomy']["NCBI_taxonomy_id"]=="83333":
            if CARD[id]["model_sequences"]["sequence"][seqid]['NCBI_taxonomy']["NCBI_taxonomy_id"] in["573","83333","287", "470"
:
                #Select E. coli K12
                k=k+1
                Oneset.add(CARD[id]["model_sequences"]["sequence"][seqid]['dna_sequence']['accession'])



        #Oneset.add(CARD[id]["model_sequences"]["dna_sequence"]["accession"])
    except:
        print id
print Oneset
print "Genes in E. coli form CART:", k
print "Genes in CART:", i

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A constant-space parser for the GeneOntology OBO v1.2 format

Version 1.0
"""
from collections import defaultdict

__author__    = "Uli Koehler"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__   = "Apache v2.0"

def processGOTerm(goTerm):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(goTerm) #Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.iteritems():
        if len(value) == 1:
            ret[key] = value[0]
    return ret

def parseGOOBO(filename):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each
    Keyword arguments:
        filename: The filename to read
    """
    with open(filename, "r") as infile:
        currentGOTerm = None
        for line in infile:
            line = line.strip()
            if not line: continue #Skip empty
            if line == "[Term]":
                if currentGOTerm: yield processGOTerm(currentGOTerm)
                currentGOTerm = defaultdict(list)
            elif line == "[Typedef]":
                #Skip [Typedef sections]
                currentGOTerm = None
            else: #Not [Term]
                #Only process if we're inside a [Term] environment
                if currentGOTerm is None: continue
                key, sep, val = line.partition(":")
                currentGOTerm[key].append(val.strip())
        #Add last term
        if currentGOTerm is not None:
            yield processGOTerm(currentGOTerm)
out = open("Antibitic_resitant_genesONLY_ARO_terms_antibiotics.txt", "wb")
out.write("ARO TERM\tAntibitics\tMechanism\n")
for goTerm in parseGOOBO("aro.obo"):
    if "relationship" in goTerm:
        #print type(goTerm["relationship"])
        if type(goTerm["relationship"])is not list:
            relationships = [goTerm["relationship"]]
        else:
            relationships = goTerm["relationship"]
        #print relationships
        for relationship in relationships:
            if re.findall("confers_resistance_to",relationship)>-1 or re.findall("confers_resistance_to_drug",relationship)>-1: #or re.findall("targeted_by",relationship)>-1 or re.findall("targeted_by_drug",relationship)>-1:
                antibiotic = re.findall("! (.*)$", relationship)
                out.write(goTerm["id"]+"\t"+antibiotic[0]+"\t"+goTerm["def"]+"\n")
out.close()