library(tidyverse)

table <- read.table("~/Bioinformatics/TF_KO/TR1/data/Table1.txt")$V1

p0 <-c("Ccr5, Cxcr3, Icos, Tnfrsf18, Havcr2, Lag3, Ifng, Il10, Ahr, Id2, Lilrb4a, Maf, Prdm1, Gzmb, Entpd1")

p <- gsub(" and" , ",", p0)

p <- strsplit(p, ", ")[[1]]

p.f <- factor(p, levels = table, ordered=TRUE)
sort(p.f)
writeClipboard(paste(sort(p.f), collapse = ", "))
