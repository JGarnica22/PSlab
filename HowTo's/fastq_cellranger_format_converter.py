""""
Date: 19 - 10 - 2020
Author: Joel Moro Borrego, IDIBAPS PSLab
Description: Conversion of CNS fastq file nomenclature to propper input format for 
cellranger performance.
"""


import os
import sys
import argvparse

def read_table(filename, fl):
    """
    Arguments: file name and fastqlist (fl)
    Description: Reads the table containing the references and 
    compares with a list of fastq identifiers to obtain  a dictionary
    mapping the conversion with the original id
    Output: Conversion dictionary
    """

    with open(filename) as f:
        conversor = {}
        for line in f:
            line = line.strip().split()
            if line and line[-1] in fl:
                g_name = line[8].replace("_", "")
                conversor[line[-1]+"_1.fastq.gz"] = g_name+"_S1_L00"+line[1]+"_R1_001.fastq.gz"
                conversor[line[-1]+"_2.fastq.gz"] = g_name+"_S1_L00"+line[1]+"_R2_001.fastq.gz"
        return conversor

def obtain_fastq_names(dir):
    """
    Arguments: directory containing the fastqfiles
    Description: Obtains the identifier of the fasrqfile
    Output: Returns a list of fastq identifiers
    """
    fastq_list = []
    for fn in os.listdir(dir):
        if fn.endswith("fastq.gz"):
            fastq_list.append(fn[:-11])

    return fastq_list   



#Set the complete path to your table as well as the file name:
table_file = sys.argv[1]

#Set the complete path to your fastq files
fastq_path = sys.argv[2]

fl = obtain_fastq_names(fastq_path)
conv = read_table(table_file, fl)


for fqn in conv:
    os.rename(fastq_path+fqn, fastq_path+conv[fqn])


# Set up -h

parser=argparse.ArgumentParser(
    description='''Conversion of CNS fastq file nomenclature to propper input format for 
cellranger performance. ''',
    epilog="""All is well that ends well.""")
parser.add_argument('-t', type=int, help="Full path to table")
parser.add_argument('-f', nargs='*', help="Full path to fastq dir")
args=parser.parse_args()
