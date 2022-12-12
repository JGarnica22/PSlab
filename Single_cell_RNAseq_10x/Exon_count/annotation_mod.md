---
title: "Annotation_file_modification"
author: "Josep Garnica"
date: "2022-12-09"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_knit$set(root.dir = '~/Bioinformatics/Cellranger_variants')
knitr::opts_chunk$set(echo = TRUE)

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
```
# Modify annotation GTF file
## load file
We talke the annotation file from the cellranger reference refdata-gex-mm10-2020-A, so files should be filtered and ready for `cellranger count`. 
```{r}
gtf <- import("genes.gtf")
gtf_df <- gtf %>% as.data.frame()
```
## Floxed genes
```{r}
# list of genes to be floxed
li <- c("Irf4", "Maf", "Bcl6", "Ahr","Cebpa",
        "Pax5", "Tcf7", "Tbx21", "Prdm1", "Pou2af1",
        "Nfia", "Tox2", "Gata3", "Id3", "Stat4",
        "Batf", "Stat3", "Atf6", "Bhlhe40", "Egr2",
        "Elk4")
```

## Loop to convert exons to genes
```{r}
new <- list()
for (i in li){
  print(paste("Converting file for entries for", i))
  
  idx <- gtf_df$gene_name == i
  df <- gtf_df[which(idx),] # obtain entries of the gene
  gtf_df <- gtf_df[-which(idx),] # remove entries from original df
  
  # change exon to gene in subset
  df_g <- df %>% dplyr::filter(type == "gene")
  df_t <- df %>% dplyr::filter(type == "transcript" &
                                 tag=="CCDS")
  df <- df %>% dplyr::filter(tag == "CCDS" &
                             !is.na(exon_number) &
                             type == "exon") #obtain only consensus exons
 
   print(paste0(length(unique(df$transcript_name)),
                " CCDS transcripts versions found for ", i,
                " : ",
                paste(unique(df$transcript_name), collapse = ", ")))
  df$type <- "gene"
  df <- df %>% mutate(gene_id=paste0(gene_id,".",exon_number),
                      gene_name=paste0(gene_name,".",exon_number),
                      mgi_id=paste0(mgi_id,".",exon_number),
                      havana_gene=paste0(havana_gene,".",exon_number
                      )
                      ) %>% 
              # in case more than exon are in CCDS keep the one wider
              distinct(width, gene_id, .keep_all = T) %>%
              arrange(gene_id, desc(width))
  # repeated exon entries to be changed for "exon" for type
  df[duplicated(df$gene_id),"type"] <- "transcript"
  df[df$type=="gene",which(is.na(df_g))] <- NA
  df[df$type=="transcript",which(is.na(df_t[1,]))] <- NA
  print(paste0(nrow(df), " exons found for ", i))
  
  # append modified entries to whole df
  gtf_df <- rbind(gtf_df, df)
  new[[i]] <- df
}

writexl::write_xlsx(new, "New_entries_annotation.xlsx")

```




# Cre recombinase gene
NC_005856.1:436-1467 Enterobacteria phage P1, complete genome

## Add Cre sequence to genome
```{bash}
cat cre.fa >> genome.fa

```

Count how many bases there are in transgene to indicate in GTF file
```{bash}
cat Cre.fa | grep -v "^>" | tr -d "\n" | wc -c

```

## Add Cre annotation to GTF
```{r}
cre_df <- data.frame(matrix(ncol=ncol(gtf_df)))
names(cre_df) <- names(gtf_df)
cre_df[1,] <- c("Cre", 1, 1032, ".", "+", ".", "gene", NA, NA,
                "NC_005856.1:436-1467", NA, "protein_coding", "Cre",
                rep(NA, ncol(gtf_df)-13))

gtf_df <- rbind(gtf_df, cre_df)
```

### Export GTF file
```{r}
gtf2 <- makeGRangesFromDataFrame(gtf_df,
                                 keep.extra.columns = T)
export(gtf2, "genes_Cre.gtf")
```

Check out exported format and content is ok.
```{bash}
grep "Irf4\>" genes_Cre.gtf
tail genes_Cre.gtf

```

Then in the cluster run the following code using `cellranger mkref`to create a new genome reference.
```{bash}
module load cellranger/6.1.2

cellranger mkref \
  --genome=custom_refdata-gex-mm10-2020-A \
  --fasta=files/genome_Cre.fa \
  --genes=files/genes_Cre.gtf

```

Once done, proceed to run `cellranger count` with the new refernece as always.
