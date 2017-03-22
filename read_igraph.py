import igraph
inp = open("Antibiotics_OZ_WT_norm_log2FC.txt")
header = inp.readline()
header = header.strip()
header = header.split("\t")
print len(header)
out= open("Antibiotics_OZ_WT_norm_log2FC.ncol","wb")
for line in inp:
    line=line.strip()
    line=line.split("\t")
    k = 1
    print len(line)
    while k <len(line):
        out.write(line[0]+"\t"+ header[k]+"\t"+line[k]+"\n")
        k = k+1
out.close()