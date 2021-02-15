# This script is to generate Hic-Pro annotation files:
# Bed file of the restriction fragments after digestion and table file of chromosomes' size.
# Please be sure that the chromosome names are the same than the ones used in your bowtie indexes!
#It was created for R version 4.0.3 (2020-10-10)

bioc.packages <- c("BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Mmusculus.UCSC.mm10", "HiTC", "rtracklayer")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
  library(i, character.only = T)
}

setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/HICPRO")

# get chr sizes, let's do it with and without confirmed sequences

chr_all <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
chr <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)[1:22] #get only confirmed chromosomes sequences, not random or others

for (i in 1:length(grep("chr", names(.GlobalEnv), value=TRUE))){
  
size <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[eval(as.symbol(grep("chr", names(.GlobalEnv), value=TRUE)[i]))]
write.table(size, file=paste0("chrom_mm10", grep("chr", names(.GlobalEnv), value=TRUE)[i] ,".sizes"), quote=FALSE, col.names=FALSE, sep="\t")

# Get  restriction fragments
# Another way to generate the list of restriction fragments, apart from digest_genome.py utility, is to use the HiTC BioConductor package.
# Note that this method only works for the genomes which are already available in BioConductor and for one restriction enzyme.

# In this case we cut for "GATC" pattern that corresponds to DpnII enzyme and
# to mm10 genome version from UCSC. Change resSite and genomePack to change enzyme and/or genome to use.

resFrag <- getRestrictionFragmentsPerChromosome(resSite="GATC", chromosomes=eval(as.symbol(grep("chr", names(.GlobalEnv), value=TRUE)[i])), 
                                                overhangs5=1, genomePack="BSgenome.Mmusculus.UCSC.mm10") 
allRF <- do.call("c",resFrag)
names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
export(allRF, format="bed", con=paste0("DpnII_resfrag_mm10", grep("chr", names(.GlobalEnv), value=TRUE)[i], ".bed"))
}

