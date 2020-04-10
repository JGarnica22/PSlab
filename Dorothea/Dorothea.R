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
# Indicate name of txt file containing expression data (raw counts)
expfile <- "Partek_TFH_Raw_counts.txt"
#Indicate name of txt file with counts normalized for sequencing depth and transformed to log2 scale
rlog_norm <- "rlog_normalized.txt"
#Indicate differential expressionl data:
DEmatrix <- "DE_Tetvs.Th0"
# Indicate populations of interest, to be compared. First indicate control population, then sample(s):
pop <- c("Th0", "TFH", "Tet")
#Indicate confidence to use:
confi2 <- c("A","B","C")
#Load or indicate your regulon as rdata or as csv
load('data/TOP10score_viperRegulon.rdata')
regucsv <- "dorothea_regulon_mouse_v1.csv"
#if regulon is not available in viper format, transform from csv
if (exists("viper_regulon_mouse")==F) {
  regulon_mouse <- read.csv(paste0("data/",regucsv), header = T, sep = ",", quote = "")
    # This function converts a regulon in data frame format to the viper format:
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
  viper_regulon_mouse <- df2regulon(regulon_mouse)
  }  
# Indicate "human" or "mouse" experiment
species <- "mouse"
# Read raw counts:
counts <- read.table(paste0("data/",expfile),
                     header = T,
                     sep = "\t",
                     quote = "",
                     dec = ".")
#Read rlog normalized counts
rlog.norm.counts <- read.table(paste0("data/",rlog_norm), header = T, sep = "\t", dec = ".", quote = "")
#if not available, make from raw counts:
if (exists("rlog.norm.counts")==F){
  row.names(counts) <- counts$gene_name
  # Eliminate unnecessary columns, we only want gene_name and counts for the different populations.
  counts2 <- data.frame(matrix(ncol=0,nrow = nrow(counts)))
  row.names(counts2) <- row.names(counts)
  for (p in c(1:length(pop))){
  counts1 <- counts[, c(grep(pop[p], colnames(counts)))]
  names(counts1) <- c(paste0(rep(pop[p], length(grep(pop[p],colnames(counts1)))),
                             1:length(grep(pop[p],colnames(counts)))))
  counts2 <- cbind(counts2,counts1)
  }
 cell_type <- NULL
 replicates <- NULL
 for (o in c(1:length(pop))) {
   ct <- rep(pop[o], length(grep(pop[o], colnames(counts))))
   cell_type <- append(cell_type, ct)
   rpli <- 1:length(grep(pop[o],colnames(counts)))
   replicates <- append(replicates,rpli)
 }
   sampe_info <- data.frame(cell_type, replicates, row.names = names(counts2))
  counts2_r <- apply(counts2, c(1,2), round)
  counts2_i <- apply(counts2_r, c(1,2), as.integer)
  dds <- DESeqDataSetFromMatrix (countData = counts2_i,
                                 colData = sample_info,
                                 design = ~ cell_type)
  dds <- dds[ rowSums (counts(dds)) > 10, ]
  DESeq.ds_f <- estimateSizeFactors(dds)
  sizeFactors(DESeq.ds_f)
  DESeq.rlog <- rlog(DESeq.ds_f, blind = TRUE)
  rlog.norm.counts <- assay (DESeq.rlog)
  }
# Clean TF names because some genes are present more than once as -1/2/...:
names(viper_regulon_mouse) <- sapply(strsplit(names(viper_regulon_mouse), split = ' - '), head, 1)
# Compute single-sample TF activities from a normalized gene expression matrix
# Estimate TF activities
TFact_gexpr <- viper(eset = E, regulon = viper_regulon_mouse, nes = T, method = 'none', minsize = 4, eset.filter = F)
pheatmap(TFact_gexpr, scale = "row", show_rownames = FALSE, title = "Viper regulon mouse")
# Save results
write.csv(TFact_gexpr, file = 'output/TFactivities_rlog_geneexpression.csv')
pdf(file = "figs/heatmap TF active gene expression.pdf", width = 7, height = 5)
print(pheatmap(TFact_gexpr, scale = "row", show_rownames = FALSE))
dev.off()

#PROBABLY NOT CORRECT!! ASK
#Perform linear regression to calculate enrichment of TFs:
controls <- sample_info$cell_type == pop[1] 
result <- apply(TFact_gexpr, 1, function(x) {
  broom::tidy(lm(x ~ !controls)) %>%
    filter(term == "!controlsTRUE") %>%
    select(-term)
})
ttest <- mutate(TF=names(result), bind_rows(result))
write.table(ttest, "output/DoRothEA_ttest.txt", dec = ".", sep = "\t", quote = FALSE)

#adjust pvalue if this presents a normal distribuition. Normal distribution visualization and saphirot test
ggqqplot(ttest$p.value)
if (shapiro.test(ttest$p.value)[2] > 0.05){
ttest$fdr <- p.adjust(ttest$p.value, method = "fdr")
ttest <- as.data.frame(ttest)
ttest_sig <- subset(ttest, ttest$fdr<= 0.05)
}

#Don't do now: requires objectes prepared later on !!
TF.sig.emat <- ttest_sig[,"TF"]
TF.sig.emat <- sort(TF.sig.emat)
TF.sig.DE <- as.vector(as.character(TF_sig[,"Regulon"]))
TF.sig.DE <- sort(TF.sig.DE)
length(TF.sig.emat)
length(TF.sig.DE)
table(TF.sig.DE %in% TF.sig.emat)

#################################################################################################
#In order to clean up and not use the whole TF elements in DoRothEA, I will try to use only
#TFs with high confidence, indicate confidence to use:
#Set confidence to subset:
regulon_mouse_confi2 <- data.frame(matrix(ncol= ncol(regulon_mouse), nrow=0))
names(regulon_mouse_confi2) <- names(regulon_mouse)
for (u in confi2) {
rm <- subset(regulon_mouse, regulon_mouse$confidence == u)
regulon_mouse_confi2 <- rbind(regulon_mouse_confi2,rm)
}
viper_regulon_mouse_confi <- df2regulon(regulon_mouse_confi2)
TFact_gexpr_confi = viper(eset = E, regulon = viper_regulon_mouse_ABC, nes = T, method = 'none',
                        minsize = 4, eset.filter = F)
pheatmap(TFact_gexpr_confi, scale = "row", show_rownames = FALSE)

#Introduce level of confidence to each TF:
TFact_gexpr_hiconf <- merge(x = TFact_gexpr_confi,
                 y = regulon_mouse_confi2, 
                 by.x = "row.names", by.y = "tf")
TFact_gexpr_hiconf <- unique(TFact_gexpr_hiconf)
TFact_gexpr_hiconf <- TFact_gexpr_ABC_conf[order(TFact_gexpr_ABC_conf$confidence), ]
write.table(TFact_gexpr_hiconf, file = paste0("output/DoRothEA_TFactivities_gexpr_", paste0(confi2, collapse = ""),".txt"), dec = ".", sep = "\t",
            quote = FALSE)
write_xlsx(TFact_gexpr_hiconf, paste0("output/DoRothEA_TFactivities_gexpr_", paste0(confi2, collapse = ""),".xlsx"))

# 2B.Compute TF activity changes from a differential gene expression signature
################################################################################
# Load differential expression signature
DEsignature <- read.table(paste0("data/", DEmatrix ,".txt"), header = T, sep = "\t",
                          dec = ".", quote = "", row.names = NULL)
if (exists("DEsignature")==F){
  counts <- read.table(paste0("data/", expfile),
                       header = T,
                       sep = "\t",
                       quote = "", 
                       dec = ".")
    row.names(counts) <- counts$gene_name
  counts2 <- data.frame(matrix(ncol=0,nrow = nrow(counts)))
  row.names(counts2) <- row.names(counts)
  for (p in c(1:length(pop))){
    counts1 <- counts[, c(grep(pop[p], colnames(counts)))]
    names(counts1) <- c(paste0(rep(pop[p], length(grep(pop[p],colnames(counts1)))),
                               1:length(grep(pop[p],colnames(counts)))))
    counts2 <- cbind(counts2,counts1)
  }
  cell_type <- NULL
  replicates <- NULL
  for (o in c(1:length(pop))) {
    ct <- rep(pop[o], length(grep(pop[o], colnames(counts))))
    cell_type <- append(cell_type, ct)
    rpli <- 1:length(grep(pop[o],colnames(counts)))
    replicates <- append(replicates,rpli)
  }
  sample_info <- data.frame(cell_type, replicates, row.names = names(counts2))
  counts2_r <- apply(counts2, c(1,2), round)
  counts2_i <- apply(counts2_r, c(1,2), as.integer)
    dds <- DESeqDataSetFromMatrix (countData = counts2_i,
                                 colData = sample_info,
                                 design = ~ cell_type)
    dds <- dds[ rowSums (counts(dds)) > 10, ]
  dds$cell_type <- relevel(dds$cell_type, ref = pop[1])
  dds <- DESeq(dds)
  for (t in c(2:length(pop))){
    resLFC <- lfcShrink(dds, coef=resultsNames(dds)[t], type="apeglm")
    resOrdered <- resLFC[order(resLFC$log2FoldChange, decreasing = T),]
    write.table(resOrdered, file = paste0("data/DE_",strsplit(resultsNames(dds)[t], "_")[[1]][3], "vs.",
                                          strsplit(resultsNames(dds)[t], "_")[[1]][5],".txt"),
                                          quote = FALSE, sep = "\t", dec = ".")
  }
}
DEsignature <- read.table(paste0("data/", DEmatrix ,".txt"), header = T, sep = "\t",
                          dec = ".", quote = "", row.names = NULL)
DEsignature$Gene <- DEsignature[,1]
DEsignature <- DEsignature[,c(7,2:6)]

# Exclude probes with unknown or duplicated gene symbol
DEsignature = subset(DEsignature, Gene != "" )
DEsignature = subset(DEsignature, ! duplicated(Gene))
# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics = matrix(DEsignature$log2FoldChange, dimnames = list(DEsignature$Gene, 'logFC') )
myPvalue = matrix(DEsignature$pvalue, dimnames = list(DEsignature$Gene, 'pvalue') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]
# Estimate TF activities
mrs = msviper(ges = mySignature, regulon = viper_regulon_mouse, minsize = 4, ges.filter = F)
TFact_DE = data.frame(Regulon = names(mrs$es$nes),
                           Size = mrs$es$size[ names(mrs$es$nes) ], 
                           NES = mrs$es$nes, 
                           p.value = mrs$es$p.value, 
                           FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TFact_DE = TFact_DE[ order(TFact_DE$p.value), ]

TFact_DE_conf <- merge(x = TFact_DE,
                     y = regulon_mouse, 
                     by.x = "row.names", by.y = "tf")
TFact_DE_conf <- TFact_DE_conf[,c("Regulon", "NES", "p.value", "FDR", "confidence")]
TFact_DE_conf <- unique(TFact_DE_conf)
TFact_DE_conf <- TFact_DE_conf[order(TFact_DE_conf$FDR), ]

TFact_DE_conf_sig <- TFact_DE_conf[which(TFact_DE_conf$FDR<=0.05),]
TFact_DE_conf_sig <- TFact_DE_conf_sig[order(TFact_DE_conf_sig$FDR), ]
table(TFact_DE_conf_sig$confidence=="A")
table(TFact_DE_conf_sig$confidence=="B")
table(TFact_DE_conf_sig$confidence=="C")

# Save results
write.csv(TFact_DE_conf, file = paste0("output/DoRothEA_TFactivities_",DEmatrix,"_allconf.csv"))
write.csv(TFact_DE_conf_sig, file = paste0("output/DoRothEA_TFactivities_",DEmatrix,"_allconf005SIG.csv"))

#Volcano NES vs FDR on TFact_DE_conf
Vol <- data.frame(NES= TFact_DE_conf$NES,
                  FDR= -log10(TFact_DE_conf$FDR), Confidence=TFact_DE_conf$confidence,
                  row.names = rownames(TFact_DE_conf))
# Vol$S <- 0
# Vol[which(Vol$FDR<2), "S"] <- "NS"
# Vol[which(Vol$NES>2 & Vol$S!="NS"), "S"] <- "UP"
# Vol[which(Vol$NES<(-2) & Vol$S!="NS"), "S"] <- "DOWN"
# Vol[which(Vol$S==0), "S"] <- "NS"
# Vol$S <- as.factor(Vol$S)
# Vol[which(Vol$FDR>=99.5), "FDR"] <- 99.5
pdf(file = paste0("figs/Volcano_", "TFact_DE_conf", ".pdf"), width = 4, height = 5)
print(ggplot(Vol, aes(x=NES, y=FDR, color = Confidence)) +
        geom_point (aes(color = factor(Confidence)), size=10, alpha = 1/3,
                    show.legend = T) +
        geom_vline(xintercept = c(-2, 2), linetype = 1, size = 0.3, col = "grey20") +
        geom_hline(yintercept = 2, linetype = 1, size = 0.3, col = "grey20") +
        scale_color_manual(values = c("red", "orange", "yellow","green", "blue")) +
        theme_light() +
        ggtitle(paste0(DEmatrix)) +
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
#boxplot NES vs FDR on TFact_DE_conf
box <- data.frame(FDR= -log10(TFact_DE_conf$FDR), Confidence=TFact_DE_conf$confidence,
                  row.names = rownames(TFact_DE_conf))
ggplot(box, aes(x=Confidence, y=FDR)) + geom_point(color = "navyblue", position = position_jitter(w = 0.3, h = 0.0)) +
  xlab("Confidence") + ylab("-log10FDR") + labs(subtitle="Line set at FDR=0.001")+
  ggtitle(paste0(DEmatrix)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face="bold"), axis.title.x = element_text(size=14, face= "bold"),
        axis.title.y = element_text(size=14, face= "bold"), plot.subtitle = ) +
  theme(axis.text.x = element_text(angle = 0, size = 14, hjust =1, face="bold"))+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "grey",
                                        size = 0.3, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey"))+
  geom_hline(yintercept = -log10(0.001), linetype = 1, size = 0.3, col = "grey20")



############################  INCLUDE THIS CODE BELOW???  ########################################


#We get too many significant TFs, I will use only TFs with confidence A to do the analysis:
regulon_mouse_A <- regulon_mouse[which(regulon_mouse$confidence=="A"),]
viper_regulon_mouse_A <- df2regulon(regulon_mouse_A)
names(viper_regulon_mouse_A) = sapply(strsplit(names(viper_regulon_mouse_A), split = ' - '), head, 1)
mrs_A = msviper(ges = mySignature, regulon = viper_regulon_mouse_A, minsize = 4, ges.filter = F)
TFact_DE_A = data.frame(Regulon = names(mrs_A$es$nes),
                             Size = mrs_A$es$size[ names(mrs_A$es$nes) ],
                             NES = mrs_A$es$nes,
                             p.value = mrs_A$es$p.value,
                             FDR = p.adjust(mrs_A$es$p.value, method = 'fdr'))
TFact_DE_A = TFact_DE_A[ order(TFact_DE_A$p.value), ]
TFact_DE_A

TFact_DE_A_sig <- TFact_DE_A[which(TFact_DE_A$FDR<=0.05),]
TFact_DE_A_sig
#From TFs with confidence A, I can only obtain 1 significant TF.

#I will use TFs with confidence A, B or C to do the analysis:
regulon_mouse_ABC <- regulon_mouse[which(regulon_mouse$confidence=="A" |
                                           regulon_mouse$confidence=="B" |
                                           regulon_mouse$confidence=="C"),]
viper_regulon_mouse_ABC <- df2regulon(regulon_mouse_ABC)
mrs_ABC = msviper(ges = mySignature, regulon = viper_regulon_mouse_ABC, minsize = 4, ges.filter = F)
TFact_DE_ABC = data.frame(Regulon = names(mrs_ABC$es$nes),
                             Size = mrs_ABC$es$size[ names(mrs_ABC$es$nes) ],
                             NES = mrs_ABC$es$nes,
                             p.value = mrs_ABC$es$p.value,
                             FDR = p.adjust(mrs_ABC$es$p.value, method = 'fdr'))
TFact_DE_ABC = TFact_DE_ABC[ order(TFact_DE_ABC$p.value), ]
TFact_DE_ABC

TFact_DE_ABC_conf <- merge(x = TFact_DE_ABC,
                       y = regulon_mouse, 
                       by.x = "row.names", by.y = "tf")
TFact_DE_ABC_conf <- TFact_DE_ABC_conf[,c("Regulon", "NES", "p.value", "FDR", "confidence")]
TFact_DE_ABC_conf <- unique(TFact_DE_ABC_conf)
TFact_DE_ABC_conf = TFact_DE_ABC_conf[ order(TFact_DE_ABC_conf$p.value), ]

write.csv(TFact_DE_ABC_conf, file = 'DoRothEA_TFactivities_DE_ABC.csv')

TFact_DE_ABC_conf_sig <- TFact_DE_ABC_conf[which(TFact_DE_ABC_conf$FDR<=0.05),]
TFact_DE_ABC_conf_sig

#In this case we can again obtain only 1 significant TF with confidence A, and only 5 TFs with confidence
#C. This is suspicious as we had before more significant TFs when using all TFs (all confidences) for the 
#analysis. If we check for less comparisons (TFs) in pval adjustment (to FDR), correction by multiple 
#testing should be less, so we would expect to find more significant results, not less.
#We will try to do pval adjustment manually and compare with the correction that is done within
#the viper function to see if the problem could be the p-values adjustment within viper. 
#For that we use the `p.adjust` function.
#We first need to create a vector with all the p-values:
pval.A <- TFact_DE_A$p.value
names(pval.A) <- row.names(TFact_DE_A)
FDR.A <- p.adjust(pval.A, method = "fdr")
#Results are the same than the obtained from viper, only 1 TF is significant looking at the FDR

#We will try the same correction for the whole set of TFs:
pval.all <- TFact_DE$p.value
names(pval.all) <- row.names(TFact_DE)
length(pval.all)
FDR.all <- p.adjust(pval.all, method = "fdr")
FDR.sig <- FDR.all<=0.05
table(FDR.sig)
#Manual correction gives no different results. The problem is with the FDR correction per se, not viper.

#We will check distribution of p-values, because if the distribution is not the usual downwards
#exponential curve, the FDR correction should not be applied.
#p-values distribution of the analysis for all TFs (all confidences A-E)
hist(TFact_DE$p.value, main = "pval distribution All-conf TFs", xlab = "pvalues")
#p-values distribution of the analysis for TFs with confidence A-C:
hist(TFact_DE_ABC$p.value, main = "pval distribution ABC-conf TFs", xlab = "pvalues")
#p-values distribution of the analysis for TFs with confidence A:
hist(TFact_DE_A$p.value, main = "pval distribution A-conf TFs", xlab = "pvalues")
#The distribution of the pvalues for only TFs with confidence A is not normal. We should not
#correct for FDR.
#It is not normal even when you check how the pvalues for the set of A-conf TFs within the whole
#analysis are distributed:
hist(TFact_DE_conf[which(TFact_DE_conf$confidence=="A"),"p.value"], 
     xlab = "pvalues",
     main = "pval distribution A-conf TFs \n within All-conf list")

#In order to check that it is only a matter of pvalue correction, we need to check also that the
#NES and the pvalues obtained for A-TFs is the same when we do it for all TFs or only for the 
#A set (in theory, viper does the analysis individually for each set):

#1. Compare NES:
NES.fil <- as.data.frame(TFact_DE_conf[which(TFact_DE_conf$confidence=="A"),"NES"])
row.names(NES.fil) <- TFact_DE_conf[which(TFact_DE_conf$confidence=="A"), "Regulon"]
names(NES.fil) <- "NES_all"

NES.A <- as.data.frame(TFact_DE_A[,"NES"])
row.names(NES.A) <- row.names(TFact_DE_A)
names(NES.A) <- "NES_A"

NES.comp <- merge(NES.A,NES.fil, by= 0)
NES.comp
#They are the exact same NES values.

#2. Compare p-values:
pval.fil <- as.data.frame(TFact_DE_conf[which(TFact_DE_conf$confidence=="A"),"p.value"])
row.names(pval.fil) <- TFact_DE_conf[which(TFact_DE_conf$confidence=="A"), "Regulon"]
names(pval.fil) <- "pval_all"

pval.A <- as.data.frame(TFact_DE_A[,"p.value"])
row.names(pval.A) <- row.names(TFact_DE_A)
names(pval.A) <- "pval_A"

pval.comp <- merge(pval.A,pval.fil, by= 0)
pval.comp
#Also the exact same p-values.
#The problem is really with FDR correction, that cannot be applied to this p-value distribution.
#In order to select interesting TFs, we will just look at the top TFs, with high NES.
