#This script is to perform DE analysis with RNAseq data
#It was created for R 3.6.2 version (2019-12-12)
#Copyright (C) 2020  Patricia Solé Sánchez
##########################################################

# Check if required packages are installed, if not install:
cran.packages <- c("ggplot2", "ggrepel")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}

bioc.packages <- c("DESeq2", "biomaRt")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

# Load packages:
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(biomaRt)


# Set your working directory (the project you are working in):
setwd("/Users/patri/Desktop/R class/Prova_GitHub_DE_analysis")

# Specify parameters to be used along the script:
## Indicate name of txt file containing expression data (raw counts)
expfile <- "Partek_TFH_Raw_counts.txt"
## Indicate populations of interest, to be compared:
#First indicate control population, then sample:
pop <- c("Th0", "TFH")


# Load data:
counts <- read.table(paste0("data/", expfile),
                     header = T,
                     sep = "\t",
                     quote = "", 
                     dec = ".")

head(counts)

#Set gene_name as rows labels
row.names(counts) <- counts$gene_name
#Eliminate unnecessary columns, we only want gene_name and counts for the different populations.
counts2 <- counts[, c(grep(pop[1], colnames(counts)), grep(pop[2], colnames(counts)))]

#Change names of columns to shorter, easier ones:
names (counts2) <- c(paste0(rep(pop[1], length(grep(pop[1],colnames(counts2)))),
                            1:length(grep(pop[1],colnames(counts2)))), 
                     paste0(rep(pop[2], length(grep(pop[2],colnames(counts2)))),
                            1:length(grep(pop[2],colnames(counts2)))))
head(counts2)

#Create conditions:
sample_info <- data.frame(cell_type = c(rep(pop[1], length(grep(pop[1], colnames(counts2)))), 
                                        rep(pop[2], length(grep(pop[2], colnames(counts2))))),
                          replicates = c(1:length(grep(pop[1],colnames(counts2))),
                                         1:length(grep(pop[2],colnames(counts2)))),
                          row.names = names (counts2))

#Transform the matrix by rounding the counts and transforming the values from "numeric" to "integer"
counts2_r <- apply(counts2, c(1,2), round)
counts2_i <- apply(counts2_r, c(1,2), as.integer)

#Generate the DESeq2 dataset:
dds <- DESeqDataSetFromMatrix (countData = counts2_i,
                               colData = sample_info,
                               design = ~ cell_type)
head(dds)
counts(dds)

#Remove genes with very low counts:
dds <- dds[ rowSums (counts(dds)) > 10, ]
counts(dds)

#DESeq’s default method to normalize read counts to account for differences in 
#sequencing depths is implemented in estimateSizeFactors:
DESeq.ds_f <- estimateSizeFactors(dds)
sizeFactors(DESeq.ds_f)


#With rlog() function from DESeq2, we directly obtain counts normalized for sequencing 
#depth and transformed to log2 scale:
DESeq.rlog <- rlog(DESeq.ds_f, blind = TRUE)
rlog.norm.counts <- assay (DESeq.rlog)
write.table(rlog.norm.counts, file = "output/rlog_normalized_DESeq.txt", quote = FALSE, sep = "\t", dec = ".")



#Plot PCA with prcomp():
#prcomp calculates principal components
pc <- prcomp(t(rlog.norm.counts))
pc_sum <- summary(pc)
PC1_varexpl <- pc_sum$importance[2,"PC1"]
PC2_varexpl <- pc_sum$importance[2,"PC2"]

pdf(file = "figs/PCA.pdf", width = 5, height = 5)
plot(x = pc$x[,1], y= pc$x[,2],
     col = c("red", "blue")[colData(DESeq.ds_f)[,1]],
     pch = 16,
     xlab = paste0("PC1 = ", round(PC1_varexpl*100,1), "%"),
     ylab = paste0("PC2 = ", round(PC2_varexpl*100,1), "%"),
     las = 1,
     main = "PCA of seq. depth normalized \n and rlog - transformed read counts ")
dev.off()

#Apply DESeq to do differential expression
dds <- DESeq(dds)
#Show dataframe from differential expression analysis. If you want, change
#the direction of comparison
res <- results(dds, contrast=c("cell_type", pop[2], pop[1]))
res

####CHECK ORDER OF COMPARISON
#Clean up genes with low counts and high variability with Shrinkage
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
head(resLFC)

#Order genes based on their pvalue
resOrdered <- resLFC[order(resLFC$log2FoldChange, decreasing = T),]


#VOLCANO PLOT:
Vol <- data.frame(FC= resOrdered$log2FoldChange,
                  sig= -log10(resOrdered$padj))
head(Vol)

plot <- ggplot(Vol, aes(x=FC, y=sig))

pdf(file = "figs/VOLCANO.pdf", width = 5, height = 5 )
plot + geom_point( aes(color=sig>2 & abs(FC)>2),
                   size = 0.3,
                   show.legend = F) + 
  scale_color_manual(values = c("TRUE" = "red","FALSE"="black")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed") +
  #if needed:
  #scale_y_continuous(breaks = c(0, 5, 10, 50, 100, 200, 250), 
  #limits =  c(0, 250), trans = "sqrt") +
  #xlim(-10, 10) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.9) +
  ggtitle(paste0(strsplit(resultsNames(dds)[2], "_")[[1]][3], " vs. ",
                 strsplit(resultsNames(dds)[2], "_")[[1]][5])) +
  xlab("log2FC") + ylab("-log10(padj)")
dev.off()

#if you need to add gene names, you may use these functions
#geom_label
#geom_text
#http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software


#If you want to show gene names for a particular set of genes:
#Change txt file to the one containing a list of genes that you want to show in volcano

# !! To unblock code, select all and press CTL+SHIFT+C
# genestoshow <- "Table1.txt"
# genes <- read.table(paste0("data/", genestoshow),
#                      sep = "\t", dec = ".", header = T)
# genes <- as.vector(genes[,1])

# set.seed(5)
# plot +
# geom_point( aes(color=sig>2 & abs(FC)>1.5),
# size = 0.3,
# show.legend = F) +
# geom_text(data= subset(Vol, rownames(Vol) %in% genes),
#                 aes(label = rownames(subset(Vol, rownames(Vol) %in% genes))),
#                 size = 2,
#                 box.padding = unit(0.15, "lines"),
#                 point.padding = NA,
#                 vjust = -0.5, nudge_x = 0.2,
#                 segment.size = 0.25) +
# scale_color_manual(values = c("TRUE" = "red","FALSE"="grey")) +
# geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", col = "grey25") +
# geom_hline(yintercept = 2, linetype = "dashed", col = "grey25") +
# theme_bw() +
# theme(panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       aspect.ratio = 0.9) +
# ggtitle(paste0(strsplit(resultsNames(dds)[2], "_")[[1]][3], " vs. ",
#                strsplit(resultsNames(dds)[2], "_")[[1]][5])) +
# xlab("log2FC") + ylab("-log10(padj)")


# BIOMART:
#select mart:
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
#Select a Dataset:
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)

#if you need to look per specific attribute:
#searchAttributes(mart = ensembl, pattern = "ensembl.*id")

#Find in Biomart ensembleID, gene type and mouse gene name
BM <- getBM(attributes=c("ensembl_gene_id","gene_biotype",
                          "external_gene_name"), mart = ensembl, verbose = T)

#Now merge this data frame with the DESeq2 results
DesBM <- merge(x= as.data.frame(resOrdered), y=BM, 
               by.x="row.names", by.y = "external_gene_name", 
               all.x = T)
write.table(DesBM, paste0("output/DESeq2_", resultsNames(dds)[2], ".txt"), 
            quote = F, dec = ".", sep = "\t")

