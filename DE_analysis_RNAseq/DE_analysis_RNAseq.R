#This script is to perform DE analysis with RNAseq data
#It was created for R 3.6.2 version (2019-12-12)
#Copyright (C) 2020  Patricia Solé Sánchez
##########################################################

# Check if required packages are installed, if not install:
cran.packages <- c("ggplot2", "ggrepel", "writexl", "pheatmap")
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
library(writexl)
library(pheatmap)
library(RColorBrewer)

# Set your working directory (the project you are working in):
setwd("/Users/patri/Desktop/R class/Prova_GitHub_DE_analysis")

# Specify parameters to be used along the script:
## Indicate name of txt file containing expression data (raw counts)
expfile <- "Partek_TFH_Raw_counts.txt"
## Indicate populations of interest, to be compared (use same name
## found in column names without numbers):
# Ex: TFH1_25887 will be TFH
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
#Eliminate unnecessary columns, we only want gene_name and counts for the different populations
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
rlog.norm.counts2 <- as.data.frame(assay (DESeq.rlog))
rlog.norm.counts2$Gene <- rownames(rlog.norm.counts2)
rlog.norm.counts2 <- rlog.norm.counts2[,c(length(rlog.norm.counts2), 1:(length(rlog.norm.counts2)-1))]
write.table(rlog.norm.counts2, "output/rlog_normalized_DESeq.txt", quote = F,
            sep = "\t", dec = ".", row.names = F)
write_xlsx(rlog.norm.counts2, "output/rlog_normalized_DESeq.xlsx")

#Plot PCA with prcomp():
#prcomp calculates principal components
pc <- prcomp(t(rlog.norm.counts))
pc_sum <- summary(pc)
PC1_varexpl <- pc_sum$importance[2,"PC1"]
PC2_varexpl <- pc_sum$importance[2,"PC2"]

plot <- as.data.frame(pc$x[,1:2])
plot$sample_type <- as.character(sample_info$cell_type)

pdf(file = paste0("figs/PCA_", pop[2], "_vs_", pop[1], ".pdf"), width = 5.5, height = 5)
ggplot(plot, aes(x = PC1, y= PC2, 
                 color = sample_type,
                 fill= sample_type)) +
  geom_point(size= 5,
             pch = 21) +
  scale_color_manual(values = c("royalblue3","green3")) + #border color, change colors if you like
  scale_fill_manual(values = c("steelblue2","greenyellow")) + #fill color, change colors if you like
  # theme_minimal() +
  #or
  theme_light() +
  labs(title= "PCA of seq. depth normalized \nand rlog - transformed read counts",
       x= paste0("PC1 (", round(PC1_varexpl*100,1), "%)"),
       y= paste0("PC2 (", round(PC2_varexpl*100,1), "%)")) +
  theme(legend.title = element_text(size=10, 
                                    face="bold")) +
  theme(#legend.background = element_rect(linetype="solid", 
    #                                 colour ="grey"),
    legend.title=element_blank())
dev.off()

#Set control population for the DE comparison
dds$cell_type <- relevel(dds$cell_type, ref = pop[1])
#Apply DESeq to do differential expression
dds <- DESeq(dds)
#Show dataframe from differential expression analysis
res <- results(dds)

#Clean up genes with low counts and high variability with Shrinkage
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
head(resLFC)
#Order genes based on their pvalue
resOrdered <- resLFC[order(resLFC$log2FoldChange, decreasing = T),]

# VOLCANO PLOT:
Vol <- data.frame(log2FC= resOrdered$log2FoldChange,
                  sig= -log10(resOrdered$padj), 
                  row.names = rownames(resOrdered))
Vol$S <- 0
Vol[which(Vol$sig<2), "S"] <- "NS"
Vol[which(Vol$log2FC>2 & Vol$S!="NS"), "S"] <- "UP"
Vol[which(Vol$log2FC<(-2) & Vol$S!="NS"), "S"] <- "DOWN"
Vol[which(Vol$S==0), "S"] <- "NS"
head(Vol)
Vol$S <- as.factor(Vol$S)
Vol[which(Vol$sig>=99.5), "sig"] <- 99.5

pdf(file = paste0("figs/Volcano_", pop[2], "_vs_", pop[1], ".pdf"), width = 4, height = 5)
ggplot(Vol, aes(x=log2FC, y=sig, color = S)) +
  geom_point (size = 0.25,
              show.legend = F) +
  
  geom_vline(xintercept = c(-2, 2), linetype = 1, size = 0.3, col = "grey20") +
  geom_hline(yintercept = 2, linetype = 1, size = 0.3, col = "grey20") +
  
  scale_color_manual(values = c("skyblue2", "grey60", "green3")) +
  theme_light() +
  ggtitle(paste0(pop[2], " vs. ", pop[1])) +
  xlab("Fold Change") + 
  ylab("-log10 (FDR)") +
  scale_y_continuous(breaks = c(1,20,40,60,80,100),
                     limits = c(-5, 100), expand = expansion(mult=0, add= c(1,0))) +
  scale_x_continuous(breaks = c(seq(-20,20, by = 1)),
                     labels = c(2^abs(-20:20)*sign(-20:20)),
                     guide = guide_axis(check.overlap = T)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey20", linetype=1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

#If you need to add gene names, you may use these functions
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
# ggplot(Vol, aes(x=log2FC, y=sig, color = S)) +
#   geom_point (size = 0.25,
#               show.legend = F) +
#   geom_point (data= subset(Vol, rownames(Vol) %in% genes & Vol$S != "NS"),
#               color = "black",
#               pch = 21,
#               size = 0.6) +
#   geom_vline(xintercept = c(-2, 2), linetype = 1, size = 0.3, col = "grey20") +
#   geom_hline(yintercept = 2, linetype = 1, size = 0.3, col = "grey20") +
#   geom_text_repel(data= subset(Vol, rownames(Vol) %in% genes & Vol$S != "NS"),
#                   aes(label = rownames(subset(Vol, rownames(Vol) %in% genes & Vol$S != "NS"))),
#                   color = "black",
#                   size = 2.8,
#                   vjust = -0.1) +
#   scale_color_manual(values = c("skyblue2", "grey60", "green3")) +
#   theme_light() +
#   ggtitle(paste0(pop[2], " vs. ", pop[1])) +
#   xlab("Fold Change") + 
#   ylab("-log10 (FDR)") +
#   scale_y_continuous(breaks = c(1,20,40,60,80,100),
#                      limits = c(-5, 100), expand = expansion(mult=0, add= c(1,0))) +
#   scale_x_continuous(breaks = c(seq(-20,20, by = 1)),
#                      labels = c(2^abs(-20:20)*sign(-20:20)),
#                      guide = guide_axis(check.overlap = T)) +
#   theme(axis.line = element_line(size = 0.3, colour = "grey20", linetype=1),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank())


# BIOMART:
#Select mart:
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
#Select a Dataset:
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)

#If you need to look per specific attribute:
#searchAttributes(mart = ensembl, pattern = "ensembl.*id")

#Find in Biomart ensembleID, gene type and mouse gene name
BM <- getBM(attributes=c("ensembl_gene_id","gene_biotype",
                         "external_gene_name"), mart = ensembl, verbose = T)

#Now merge this data frame with DESeq2 results
DesBM <- merge(x= as.data.frame(resOrdered), y= BM, 
               by.x= "row.names", by.y = "external_gene_name", 
               all.x = T)
DesBM <- DesBM[ , c(1, length(DesBM)-1, length(DesBM), 2:(length(DesBM)-2))]
colnames(DesBM) <- c("Gene", colnames(DesBM[2:length(DesBM)]))
DesBM <- DesBM[order(DesBM$log2FoldChange, decreasing = T),]
DesBM <- DesBM[!duplicated(DesBM$Gene),]

write.table(DesBM, paste0("output/DESeq2_", 
                          strsplit(resultsNames(dds)[2], "_")[[1]][3], "_vs_",
                          strsplit(resultsNames(dds)[2], "_")[[1]][5], ".txt"), 
            quote = F,
            sep = "\t", dec = ".", row.names = F)
write_xlsx(DesBM, paste0("output/DESeq2_", 
                         strsplit(resultsNames(dds)[2], "_")[[1]][3], "_vs_",
                         strsplit(resultsNames(dds)[2], "_")[[1]][5], ".xlsx"))

# HEATMAP:
DEgenes <- DesBM[which(abs(DesBM$log2FoldChange)>2 & DesBM$padj<0.01), "Gene"]

pdf(file = paste0("figs/Heatmap_", pop[2], "_vs_", pop[1], ".pdf"), width = 4, height = 6)
pheatmap(rlog.norm.counts[DEgenes,], scale = "row",
         rev(brewer.pal(8, ("RdBu"))), show_rownames = F,
         # cluster_rows = F, cluster_cols = F, na_col = "white",
         # border_color = "black",
         # cellwidth = 20, cellheight = 10,
         # main = "Table 1", 
         # fontsize = 10,
         # fontsize_row = 9,
         # fontsize_col = 12, angle_col = 90,
         # width = 800, height = 4000)
)
dev.off()
