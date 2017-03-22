inp = open("targets_iSNP_network_by_TFs_and_miRNAs.txt")

outdic = {}
for lin in inp:
    line = line.strip()
    line = line.split("\t")
    effectors = line[2] .split(":")
    if len(effectors)>1:
        for effector in effectors:
            if effector:
                try:
                    outdic[line[0]+"\t"+line[1]].add(effector)
                except:
                    outdic[line[0] + "\t" + line[1]] = set()
                    outdic[line[0] + "\t" + line[1]].add(effector)
    else:
        effector = effectors
        try:
            outdic[line[0] + "\t" + line[1]].add(effector)
        except:
            outdic[line[0] + "\t" + line[1]] = set()
            outdic[line[0] + "\t" + line[1]].add(effector)
print outdic

