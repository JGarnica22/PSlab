"""
Program: get_contact_table
Description: From a directory with .log files extracted from ccp4, generates a consensus
table given that the directory format is allele --> pdb --> bond_comparisons --> comparison

"""

import os

#Assign the protein structure directory:
wd = "PROTEIN_STRUCTURES/"


def read_info(fn):
    """
    Arguments: Info file from protein structures
    Description: Returns a dictionary that allows us to understand which is the comparison
    performed in each bond subdirectory logfile
    """
    conversion = {}
    with open(fn) as f:
        i = 0
        for line in f:
            line = line.split()
            if i == 2:
                for conv in line:
                    conv = conv.split("=")
                    conversion[conv[0]] = conv[1]
            i += 1
    return conversion

def check_gap(log_pos, seq, diff, ref, res):
    """
    """
    subseq = ""
    log_pos = int(log_pos)

    if log_pos < ref:
        return log_pos
    else:
        if seq[log_pos+diff] == res:
            return log_pos+diff
        elif res == "Q" and log_pos == 65:
            return log_pos+diff-1
        else:
            subseq = seq[(56+diff):(log_pos+diff)]
            c = subseq.count("-")
            if seq[log_pos+diff+c] == "-":
                temp = log_pos+diff+c
                while seq[temp] == "-":
                    temp+=1
                return temp
            else:
                if seq[log_pos+diff+c] == res:
                    return log_pos+31+c
                else:
                    return log_pos+31+c-1

def calculate_position(log_pos, res, mhc, seqsa, seqsb, chain, al, pdb, aminoacid_conv):
    """
    Arguments: All position information + aminoacid traduction dictionary
    Description: Calculates the new position according to the alginment and a reference sequence numeration.
    **IMPORTANT**: NUMBERS HAVE BEEN OBTAINED BY LOOKING MANUALLY TO PDB MHC SEQUENCES IN PYMOL, CALCULATIONS
    ARE SPECIFIC TO THE PROJECT
    """
    new_position = ""
    #Define diffs
    corr_pos = 0
    if res == "Wat":
        new_position = log_pos
    else:
        if chain == "a" and int(log_pos)>20:
            #Different for each pdb  (look in pymol xd)
            if al == "I-Ag7" and "6" in pdb:
                corr_pos = int(log_pos)+24
            elif al == "I-Ab" and "3" in pdb:
                corr_pos = int(log_pos)+25
            else:
                corr_pos = int(log_pos)+26

            if seqsa[al][corr_pos] != aminoacid_conv[res]:
                print("For allele", al, "in chain", chain, "the position", aminoacid_conv[res], log_pos, mhc, "is different from alignment position",seqsa[al][corr_pos], corr_pos)
                new_position = log_pos
            else:
                new_position = corr_pos
            
        #Different for each pdb (look in pymol xd)
        if chain == "b" and int(log_pos) > 20:

            if al == "I-Ag7" and "6D" in pdb:
                if int(log_pos) <= 65:
                    corr_pos = int(log_pos)+30
                elif int(log_pos) == 66:
                    corr_pos = 98
                else:
                    corr_pos = int(log_pos)+32
            elif al == "I-Ag7" and ("3" in pdb  or "ny" in pdb):
                if (int(log_pos) != 66 and "3" in pdb) or ("ny" in pdb and int(log_pos)<65):
                    corr_pos = int(log_pos)+31
                elif int(log_pos) == 66 or ("ny" in pdb and int(log_pos) == 65):
                    corr_pos = 98
                if "ny" in pdb and int(log_pos) > 66:
                    corr_pos = int(log_pos)+33
            elif al == "DQ8" and pdb == "6xcp":
                corr_pos = int(log_pos)+4
            elif al == "DQ8" and pdb == "6xc9":
                corr_pos = int(log_pos)+17
            elif al == "DQ8" and pdb == "6xco":
                corr_pos = int(log_pos)-11
            elif (al == "I-Ab" and "3c" in pdb):
                corr_pos = int(log_pos)+5
            else:
                corr_pos = int(log_pos)+31
            if  seqsb[al][corr_pos] != aminoacid_conv[res]:
                print("For allele", al, "in chain", chain, "the position", aminoacid_conv[res], log_pos, "is different from aligment position",seqsb[al][corr_pos], corr_pos)
                new_position = log_pos
            else:
                new_position = corr_pos-31
        else:
            new_position = log_pos
    return str(new_position)

def read_log(log, al, prot, bond, translator, seqs_a, seqs_b,):
    """
    Arguments: Log file, allele, pdb id, bond under study, translator dictionary
    Description: Returns a nested list in which each inner list corresponds to a modified line
    that will can be used for the merged list
    """
    #Define an amino converter
    aminoacid_conv = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M', "Wat": "Wat"}
    nested_list = []
    with open(log) as f:
        x = False
        for line in f:
            line = line.split()
            #get annotation
            if line and prot in line[0]:
                temp = line[0].split("_")
                chain = translator[temp[1]][0]
                chain_og = chain
                source_chain = translator[temp[1]]
                source_chain_og = source_chain
                target_chain = translator[temp[2]]
            #Start reading table
            if x and line:
                #Break when no more table to read in log file
                if line and line[0] == '<B><FONT':
                    x = False
                    break
                #Reading table for entries with source information 
                elif len(line) == 10 or len(line) == 9:
                    source_res = line[0]
                    #Last number corresponds to position difference log positions and aligned string
                    #Calculation was made apart xd
                    #source_position = line[1][:-1]
                    source_position = calculate_position(line[1][:-1], source_res, source_chain, seqs_a, seqs_b, chain_og, al, prot, aminoacid_conv)                
                    if chain == "a":
                         alter_position = str(int(source_position)-2)
                    else:
                        alter_position = source_position
                    source_res = line[0]+source_position
                    source_restrans = aminoacid_conv[line[0]]
                    source_atom = line[2]
                    target_res = line[4]
                    target_position = line[5][:-1]
                    target_atom = line[6]
                    distance = line[8]

                #Reading table for entries with source information given previously
                elif len(line) == 7 or len(line) == 6:
                    target_res = line[1]
                    target_position = line[2][:-1]
                    target_atom = line[3]
                    distance = line[5]
                #If chain beta has positions below 20, it corresponds to peptide
                if chain_og == "b" and int(source_position)<20:
                    if source_chain_og[-1] in "1234567890":
                        source_chain = "PEPT"+source_chain_og[-1]
                        chain = "P"
                    else:
                        source_chain = "PEPT"
                        chain = "P"
                #When that is not the case go with the originals
                else:
                    chain = chain_og
                    source_chain = source_chain_og
                nested_list.append([source_res, source_restrans, source_position, alter_position, source_atom, target_res, target_position, target_atom, distance, source_chain, target_chain, chain, bond, prot, al])
            if line and "source" == line[0]:
                x = True
            
    return nested_list

def get_sequence(file):
    """
    Arguments: Fasta file
    Description: Returns a dictionary with key: Seq ID and value: Sequence
    IMPORTANT: In particular for this case there is a deletion of the first 87 nucleotides
    since were not needed for this particular study.
    """
    Sequences = {}
    with open(file) as f:
        for line in f:
            if len(line) > 1:
                line = line.strip()
                if line.startswith(">"):
                    seq_name = line[1:]
                    if seq_name not in Sequences:
                        Sequences[seq_name[2:]] = ""
                else:
                    seq = line.replace('\n','')
                    Sequences[seq_name[2:]] += seq
    return Sequences        


######
## Main exec section
######

#Get sequences aligned for consensus reisude
alphas = "alpha_aligned.fasta"
betas = "beta_aligned.fasta"
seqs_a = get_sequence(alphas)
seqs_b = get_sequence(betas)

#Get dir alleles
maindir , alleles, filenames = next(os.walk(wd))
#Header for the macro merged file
header = ['source_residue', "source_aa", "source_position", "alter_position", 'source_atoms', 'target_residue', "target_position", "target_atoms", 'distance', 'source_chain', "target_chain", "chain", 'bond', 'pdb', 'allele']
merged_nest = [header]
for al in alleles:
    #For each dir allele we have pdb directory
    adir , pdirs, filenames = next(os.walk(wd+al+"/"))
    for prot in pdirs:
            #For each prot dir we have an info file and a bond type subdirectory
            pdir , bdirs, info = next(os.walk(wd+al+"/"+prot+"/"))
            #Mac sometimes creates file .DS_store...
            if ".DS_Store" in info:
                translator = read_info(wd+al+"/"+prot+"/"+info[1])
            else:
                translator = read_info(wd+al+"/"+prot+"/"+info[0])
            #Translator is  a dictionary that allows us to understand the comparisons
            for bond in bdirs:
                bdir , subdirs, logfiles = next(os.walk(wd+al+"/"+prot+"/"+bond+"/"))
                #Looping over logfiles
                for logfile in logfiles:
                    #Check that the log being read is .txt and not .html
                    if logfile.endswith(".log"):
                        #Read the log
                        log = wd+al+"/"+prot+"/"+bond+"/"+logfile
                        temp_nest = read_log(log, al, prot, bond, translator, seqs_a, seqs_b)
                        for l in temp_nest:
                            merged_nest.append(l)

### Writting out merged file

with open('merged_table', 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in merged_nest)


            
