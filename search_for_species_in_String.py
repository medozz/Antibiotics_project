import os
def wirte_out_species(species_IDs,string_file):
    species_out_dic = {}
    for species_ID in species_IDs:
        out= open(str(species_ID)+"string_data.txt","wb")
        species_out_dic[species_ID] = out
    inp =open(string_file)
    for line in inp:
        line= line.strip()
        listline = line.split(" ")
        for species_ID in species_IDs:
            out = species_out_dic[species_ID]
            if listline[0].find(str(species_ID))>-1 and listline[1].find(str(species_ID))>-1:
                line = line.replace((str(species_ID)+"."),"")
                out.write(line)
                out.write("\n")
    for species_ID in species_out_dic:
        species_out_dic[species_ID].close()

wirte_out_species([716541,1006551], "F:\Works\Databeses\STRING\Protein_links_STRING\protein.links.detailed.v10.txt")

def nodes_in_species(folder):
    for file_ in os.listdir(folder):
        if file_.endswith("string_data.txt"):
            inp = open(os.path.join(folder,file_))
            nodes = set()
            for line in inp:
                line=line.split(" ")
                nodes.add(line[0])
                nodes.add(line[1])
            out = open(file_.replace("string_data.txt", "string_nodes.txt"),"wb")
            out.write("\n".join(nodes))
            print file_, len(nodes)

nodes_in_species("C:\Working_projects\Antibiotics_project")