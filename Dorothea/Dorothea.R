#This script is to perform DE analysis with RNAseq data
#It was created for R 3.6.3 version (2020-04-04)
#Copyright (C) 2020  Patricia Sole Sanchez

# Check if required packages are installed, if not install:
cran.packages <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "car", "ggpubr" )
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}
bioc.packages <- c("DESeq2", "biomaRt", "viper")
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
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(biomaRt)
library(writexl)
library(viper)
library(magrittr)
library(purrr)
library(dplyr)
library(readxl)
library(car)
library(ggpubr)

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/Dorothea")

## Specify parameters to be used along the script:
# Indicate name of txt file with counts normalized for sequencing depth and transformed to log2 scale
rlog_norm <- "rlog_normalized_DESeq.txt"
# Indicate name of txt file containing DESeq2 analysis result previously performed
Desq2file <- "DESeq2_TFH_vs_TH0.txt"
# First indicate control population, then sample:
pop <- c("Th0", "TFH")
# Indicate species used in the analysis, human or mouse:
species <- "mouse"
# Indicate confidence levels to use:
confi2 <- c("A","B","C")


# Read rlog normalized counts:
rlog.norm.counts <- read.table(paste0("data/",rlog_norm), header = T, sep = "\t", dec = ".", quote = "")
# Read the DESeq2 results:
DESeq2 <- read.table (file = paste0("data/",Desq2file),
                      sep = "\t", 
                      quote = "",
                      dec = ".", 
                      header = T, row.names = "Gene")
# Load your regulon (table) as rda:
load(paste0("data/dorothea_", species, ".rda"))


# Prepare your normalized data to work with after:
rownames(rlog.norm.counts) <- rlog.norm.counts$Gene
rlog.norm.counts$Gene <- NULL

# Prepare differential expression signature
DESeq2$Gene <- rownames(DESeq2)
DESeq2 <- DESeq2[, c("Gene", "log2FoldChange", "padj")]
# Exclude probes with unknown or duplicated gene symbol
DESeq2 = subset(DESeq2, Gene != "" )
DESeq2 = subset(DESeq2, ! duplicated(Gene))

# Convert your network (table) to the viper regulon format:
dorothea_network$likelihood <- rep(1, nrow(dorothea_network))
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}
viper_regulon <- df2regulon(dorothea_network)
# Clean TF names because some genes are present more than once as -1/2/...:
names(viper_regulon) <- sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)


## A. Compute single-sample TF activities from a normalized gene expression matrix
#######################################################################################
# Estimate TF activities:
TFact_gexpr <- viper(eset = rlog.norm.counts, regulon = viper_regulon, nes = T, method = 'none', minsize = 4, eset.filter = F)
pheatmap(TFact_gexpr, scale = "row", show_rownames = FALSE, 
         main = paste("Viper regulon", pop[2], "vs.", pop[1], species))

TFact_gexprdf <- as.data.frame(TFact_gexpr)
TFact_gexprdf$TF <- row.names(TFact_gexprdf)
TFact_gexprdf <- TFact_gexprdf[,c(length(TFact_gexprdf), 1:(length(TFact_gexprdf)-1))]
# Introduce level of confidence to each TF:
TFact_gexprdf2 <- merge(x = TFact_gexprdf,
                        y = dorothea_network[,c("tf","confidence")], 
                        by.x = "TF", by.y = "tf")
TFact_gexprdf2 <- unique(TFact_gexprdf2)
write_xlsx(TFact_gexprdf2, paste0("output/DoRothEA_TFactivities_gexpr_", pop[2], "_", pop[1], "_allconf.xlsx"))

pdf(file = "figs/Heatmap_TF_activity_gexpr_allconf.pdf", width = 5, height = 7)
pheatmap(TFact_gexpr, scale = "row", show_rownames = FALSE, 
         main = paste("Viper regulon", pop[2], "vs.", pop[1], species))
dev.off()
#while (!is.null(dev.list()))  dev.off()


# In order to clean up and not use the whole TF elements in DoRothEA, use only selected confidences
regulon_confi2 <- data.frame(matrix(ncol= 5, nrow=0))
names(regulon_confi2) <- names(dorothea_network)
for (u in confi2) {
  rm <- subset(dorothea_network, dorothea_network$confidence == u)
  regulon_confi2 <- rbind(regulon_confi2,rm)
}
viper_regulon_confi <- df2regulon(regulon_confi2)
TFact_gexpr_confi = viper(eset = rlog.norm.counts, regulon = viper_regulon_confi, nes = T, method = 'none',
                          minsize = 4, eset.filter = F)
pheatmap(TFact_gexpr_confi, scale = "row", show_rownames = FALSE,
         main = paste0("Viper regulon ", pop[2], " vs. ", pop[1], " ", species,
                       "\nConfidence ", paste0(confi2, collapse = "")))

pdf(file = paste0("figs/Heatmap_TF_activity_gexpr_Confidence_", paste0(confi2, collapse = ""), ".pdf"),
    width = 5, height = 7)
pheatmap(TFact_gexpr_confi, scale = "row", show_rownames = FALSE,
         main = paste0("Viper regulon ", pop[2], " vs. ", pop[1], " ", species,
                       "\nConfidence ", paste0(confi2, collapse = "")))
dev.off()

# Introduce level of confidence to each TF:
TFact_gexpr_confidf <- as.data.frame(TFact_gexpr_confi)
TFact_gexpr_confidf$TF <- row.names(TFact_gexpr_confidf)
TFact_gexpr_confidf <- TFact_gexpr_confidf[,c(length(TFact_gexpr_confidf), 1:(length(TFact_gexpr_confidf)-1))]
TFact_gexpr_hiconf <- merge(x = TFact_gexpr_confidf,
                            y = regulon_confi2[,c("tf","confidence")], 
                            by.x = "TF", by.y = "tf")
TFact_gexpr_hiconf <- unique(TFact_gexpr_hiconf)
TFact_gexpr_hiconf <- TFact_gexpr_hiconf[order(TFact_gexpr_hiconf$confidence), ]
write.table(TFact_gexpr_hiconf, file = paste0("output/DoRothEA_TFactivities_gexpr_", 
                                              pop[2], "_", pop[1], "_",
                                              paste0(confi2, collapse = ""),".txt"), 
            dec = ".", sep = "\t", quote = FALSE)
write_xlsx(TFact_gexpr_hiconf, paste0("output/DoRothEA_TFactivities_gexpr_",
                                      pop[2], "_", pop[1], "_",
                                      paste0(confi2, collapse = ""),".xlsx"))


# B.Compute TF activity changes from a differential gene expression signature (GES)
#####################################################################################
# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics = matrix(DESeq2$log2FoldChange, dimnames = list(DESeq2$Gene, 'logFC') )
myPvalue = matrix(DESeq2$padj, dimnames = list(DESeq2$Gene, 'padj') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]
# Estimate TF activities
mrs = msviper(ges = mySignature, regulon = viper_regulon, minsize = 4, ges.filter = F)
TFact_DE = data.frame(Regulon = names(mrs$es$nes),
                      Size = mrs$es$size[ names(mrs$es$nes) ], 
                      NES = mrs$es$nes, 
                      p.value = mrs$es$p.value, 
                      FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TFact_DE = TFact_DE[ order(TFact_DE$p.value), ]

TFact_DE_conf <- merge(x = TFact_DE,
                       y = dorothea_network, 
                       by.x = "row.names", by.y = "tf")
TFact_DE_conf <- TFact_DE_conf[,c("Regulon", "NES", "p.value", "FDR", "confidence")]
TFact_DE_conf <- unique(TFact_DE_conf)
TFact_DE_conf <- TFact_DE_conf[order(TFact_DE_conf$FDR), ]

TFact_DE_conf_sig <- TFact_DE_conf[which(TFact_DE_conf$FDR<=0.05),]
TFact_DE_conf_sig <- TFact_DE_conf_sig[order(TFact_DE_conf_sig$FDR), ]

# Save results
write_xlsx(TFact_DE_conf, paste0("output/DoRothEA_TFactivities_GES_", pop[2], "_vs_", pop[1],"_allconf.xlsx"))
write_xlsx(TFact_DE_conf_sig, paste0("output/DoRothEA_TFactivities_GES_", pop[2], "_vs_", pop[1],"_allconf_FDR005.xlsx"))

# Volcano NES vs FDR on TFact_DE_conf
Vol <- data.frame(NES= TFact_DE_conf$NES,
                  FDR= -log10(TFact_DE_conf$FDR), Confidence=TFact_DE_conf$confidence,
                  row.names = rownames(TFact_DE_conf))

pdf(file = paste0("figs/Volcano_", "TFact_", pop[2], "_vs_", pop[1], ".pdf"), width = 6, height = 5)
print(ggplot(Vol, aes(x=NES, y=FDR, color = Confidence)) +
        geom_point (aes(color = factor(Confidence)), size=6.5, alpha = 1/3,
                    show.legend = T, position = position_jitter(w = 1, h = 3.0)) +
        geom_vline(xintercept = c(-2, 2), linetype = 1, size = 0.3, col = "grey20") +
        geom_hline(yintercept = 2, linetype = 1, size = 0.3, col = "grey20") +
        scale_color_manual(values = c("red", "orange", "yellow","green", "blue")) +
        theme_light() +
        ggtitle(paste(pop[2], "vs", pop[1])) +
        xlab("Normalized enrichment score") + 
        ylab("-log10 (FDR)") +
        # scale_y_continuous(breaks = c(1,20,40,60,80,100),
        #                    limits = c(-5, 100), expand = expansion(mult=0, add= c(1,0))) +
        # scale_x_continuous(breaks = c(seq(-20,20, by = 1)),
        #                    labels = c(2^abs(-20:20)*sign(-20:20)),
        #                    guide = guide_axis(check.overlap = T)) +
        theme(axis.line = element_line(size = 0.3, colour = "grey20", linetype=1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank()))
dev.off()

# Boxplot NES vs FDR on TFact_DE_conf
pdf(file = paste0("figs/Graphs_", "TFact_", pop[2], "_vs_", pop[1], ".pdf"), width = 4, height = 5)
box <- data.frame(FDR= -log10(TFact_DE_conf$FDR), Confidence=TFact_DE_conf$confidence,
                  row.names = rownames(TFact_DE_conf))
ggplot(box, aes(x=Confidence, y=FDR)) + geom_point(color = "navyblue", position = position_jitter(w = 0.3, h = 0.0)) +
  xlab("Confidence") + ylab("-log10FDR") + labs(subtitle="Line set at FDR=0.001")+
  ggtitle(paste0(pop[2], "_vs_", pop[1])) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face="bold"), axis.title.x = element_text(size=14, face= "bold"),
        axis.title.y = element_text(size=14, face= "bold"), plot.subtitle = ) +
  theme(axis.text.x = element_text(angle = 0, size = 14, hjust =1, face="bold"))+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "grey",
                                        size = 0.3, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey"))+
  geom_hline(yintercept = -log10(0.001), linetype = 1, size = 0.3, col = "grey20")
dev.off()
