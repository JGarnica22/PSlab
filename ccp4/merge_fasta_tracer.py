"""
Author: Joel Moro Borrego
Date: 26-11-2021
Description: Merges all the multifastas from TRACER in a directory as the same fasta file, while keeping track
of the original file in each fasta identifier as well as additional TRACER information
"""
import os 
#Path to tracer files
tracer_dir = "data/tracer"
#Path to clonotypes file
clonos = "data/tracer_summary/TCR_summary.txt"
#Path to sampletype file
stype = "data/sampletypes.txt"

def read_fasta(fn):
    """
    Argument: Fasta filename
    Output: Returns a dictionary key = FastaID and value = sequence.
    """
    sequences = {}
    id = ""
    #mark used to keep track of origin fasta filename
    mark = fn.split("/")[-1][:-11]
    with open(fn) as f:
        for line in f:
            line = line.split()
            if len(line) > 0 and ">" in line[0]:
                id = line[0]+"|"+mark
            elif line:
                sequences[id] = line[0]
    return sequences

def generate_clonotypes(fn):
    start = False
    i = 1
    clonotypes = {}
    with open(fn) as f:
        for line in f:
            line = line.strip(",").split()
            if start and line:
                l = line[-1]
                line = [x[:-1] for x in line[:-1]]
                line.append(l)
                for id in line:
                    clonotypes[id] = i
                i += 1
            if line and line[0] == "It":
                start = True
    return clonotypes

def read_stype(fn):
    sampletypes = {}
    header = False
    with open(fn) as f:
        for line in f:
            line = line.split()
            if header:
                sampletypes[line[0]] = [line[-2], line[-1]]
            else:
                header = True
    return sampletypes

dir , subdirs, filenames = next(os.walk(tracer_dir))
sampletypes = read_stype(stype)
clonotypes = generate_clonotypes(clonos)
a = 0
#Keep clonotypes just in case
with open('outs/clonotypes.txt', 'w') as of:
    for id in clonotypes:
        cg = clonotypes[id]
        print(id+"\t"+str(cg), file = of)
        a+=1

norecombs = []
recombs = []
sequences = {}
for sdir in subdirs:
    dir , subdirs, filenames = next(os.walk(tracer_dir+"/"+sdir+"/filtered_TCR_seqs"))
    if sdir+"_TCRseqs.fa" in filenames:
        seqs = read_fasta(tracer_dir+"/"+sdir+"/filtered_TCR_seqs/"+sdir+"_TCRseqs.fa")
        recombs.append(sdir)
        for s in seqs:
            sequences[s] = seqs[s]
    else:
        norecombs.append(sdir)

###
# Generate alpha and beta merged files 
###

tcr_map = {}
j = 1

#Write alpha
with open('outs/tracer_merged_fasta_alpha.fa', 'w') as of:
    for s in sequences:
        s_list = s.split("|")
        tcr = s_list[-1]
        #Map TCR
        if tcr not in tcr_map:
            tcr_map[tcr] = j
            j += 1
        if s_list[-1] in clonotypes:
            cgs = clonotypes[s_list[-1]]
        else:
            cgs = "U"
        if s_list[-1] in sampletypes:
            treat = sampletypes[s_list[-1]]
        else:
            print(s_list[-1], "algo pasa xd")
        if s_list[2] == "A":
            gene = s_list[-2].split("_")
            id = ">"+s_list[2]+"|"+str(cgs)+"|"+gene[0]+"_"+gene[-1]+"|"+treat[0]+"|TCR"+str(tcr_map[tcr])
            if len(gene) == 4:
                id = ">"+s_list[2]+"|"+str(cgs)+"|"+gene[0]+"-"+gene[1]+"_"+gene[-1]+"|"+treat[0]+"|TCR"+str(tcr_map[tcr])
            print(id, file = of)
            print(sequences[s], file = of)
            print("", file = of)
            


of.close()
#Write beta
with open('outs/tracer_merged_fasta_beta.fa', 'w') as of:
    for s in sequences:
        s_list = s.split("|")
        tcr = s_list[-1]
        #Map TCR
        if tcr not in tcr_map:
            tcr_map[tcr] = j
            j+=1
        if s_list[-1] in clonotypes:
            cgs = clonotypes[s_list[-1]]
        else:
            cgs = "U"
        if s_list[-1] in sampletypes:
            treat = sampletypes[s_list[-1]]
        else:
            print(s_list[-1], "algo pasa xd")
        if s_list[2] == "B":
            gene = s_list[-2].split("_")
            id = ">"+s_list[2]+"|"+str(cgs)+"|"+gene[0]+"_"+gene[-1]+"|"+treat[0]+"|TCR"+str(tcr_map[tcr])
            if len(gene) == 4:
                id = ">"+s_list[2]+"|"+str(cgs)+"|"+gene[0]+"-"+gene[1]+"_"+gene[-1]+"|"+treat[0]+"|TCR_"+str(tcr_map[tcr])
            print(id, file = of)
            print(sequences[s], file = of)
            print("", file = of)

of.close()

with open('outs/tcr_patri_map.txt', 'w') as of:
    for tcr in tcr_map:
        print(tcr+"\t"+str(tcr_map[tcr]), file = of)
of.close()


print("Recombinants:",len(recombs), "| Non recombinants:", len(norecombs),  "| Total of inputs:",len(recombs) + len(norecombs))
print("Number of fastas with clonogroup:",a)
print("Number of uniques:", 560 - a)