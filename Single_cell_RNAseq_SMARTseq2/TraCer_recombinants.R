# Script to generate comprehensive tables of TraCer results

library(tidyverse)
library(writexl)

setwd("/home/josep/Documents/Bioinformatics/Terminal/SmartSeq2/")

# import recombinants.txt file from TraCer output which includes all the productive and non-productive TCR for all analyzed cells.

recomb <- read.table("recombinants.txt", sep = "\t", quote = "",
                  dec = ".", header = T, na.strings = T, fill = T) %>%
          filter(productive == "True") %>% # select only productive TCR
          dplyr::select(-productive)

reco <- recomb %>% unite(TCR, c(locus, recombinant_id, reconstructed_length), sep = ".", remove = T) # unite all columns into one per TCR
res <- lapply(split(reco, reco$cell_name), function(i) c(t(i[ -1 ]))) # split data frame based on cell_name
res <- do.call(rbind, lapply(res, `length<-`, max(lengths(res)))) %>% as.data.frame.array()
res$cell_name <- rownames(res) 
A1 <- res[grep("A\\.", res$V1),] %>% dplyr::select(V1, cell_name)
A2 <- res[grep("A\\.", res$V2),] %>% dplyr::select(V2, cell_name)
B1 <- res[grep("B\\.", res$V1),] %>% dplyr::select(V1, cell_name) 
B2 <- res[grep("B\\.", res$V2),] %>% dplyr::select(V2, cell_name)
B3 <- res[grep("B\\.", res$V3),] %>% dplyr::select(V3, cell_name) %>% rename(V1 = V3)
B4 <- res[grep("B\\.", res$V4),] %>% dplyr::select(V4, cell_name) %>% rename(V2 = V4)
B1.1 <- rbind(B1,B3)
B2.1 <- rbind(B2,B4)


fin <- data.frame(matrix(nrow = nrow(res), ncol = 0))
fin$cell_name <- res$cell_name 
fin <- merge(fin, A1, by = "cell_name", all = T)
fin <- merge(fin, A2, by = "cell_name", all = T)
fin <- merge(fin, B1.1, by = "cell_name", all = T)
fin <- merge(fin, B2.1, by = "cell_name", all = T)
rownames(fin) <- fin$cell_name
names(fin) <- c("cell names", "A.1", "A.2", "B.1", "B.2")

fin <- fin %>% separate(A.1, into = c("locus.A1", "A1", "length.A1"), sep = "\\.", remove=T) %>% 
separate(A.2, into = c("locus.A2", "A2", "length.A2"), sep = "\\.", remove=T) %>% 
separate(B.1, into = c("locus.B1", "B1", "length.B1"), sep = "\\.", remove=T) %>% 
separate(B.2, into = c("locus.B2", "B2", "length.B2"), sep = "\\.", remove=T) %>% 
select(-c(locus.A1, locus.A2, locus.B1, locus.B2))

for (i in 1:nrow(fin)){
  if(is.na(fin[i, 6]) == TRUE) {
    fin[i,6:7] <- fin[i,8:9]
    fin[i,8:9] <- NA
    }
  }

write_xlsx(fin, "Reconstructed_TCR.xlsx", )
