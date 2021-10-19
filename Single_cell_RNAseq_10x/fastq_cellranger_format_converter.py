""""
Date: 19 - 10 - 2020
Author: Joel Moro Borrego, IDIBAPS PSLab
Description: Conversion of CNS fastq file nomenclature to propper input format for 
cellranger performance.
"""


import os


def read_table(filename, fl):
    """
    Arguments: file name and fastqlist (fl)
    Description: 
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
    Description:
    """
    fastq_list = []
    for fn in os.listdir(dir):
        if fn.endswith("fastq.gz"):
            fastq_list.append(fn[:-11])

    return fastq_list   



#Set the path to your table as well as the file name:
table_file = "2020_05_10xGenomics_BDC_INS13-21_SANTAMARIA_09/SANTAMARIA_09.txt"

#Set path to your fastq files
fastq_path = "2020_05_10xGenomics_BDC_INS13-21_SANTAMARIA_09/fastq_files_copia/"

fl = obtain_fastq_names(fastq_path)
conv = read_table(table_file, fl)


for fqn in conv:
    os.rename(fastq_path+fqn, fastq_path+conv[fqn])
