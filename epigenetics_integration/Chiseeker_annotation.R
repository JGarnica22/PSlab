annot  <- function(FOK){
  ann <- ChIPseeker::annotatePeak(FOK, tssRegion=c(-2000, 2000),
                         TxDb=TxDb, annoDb="org.Mm.eg.db", level = "gene",
                         flankDistance = 50000, addFlankGeneInfo = T)
fok <- as.data.frame(ann)
fok$type <- fok$annotation
fok$type[fok$type != "Distal Intergenic"] <- "gene_body"
fok$type[fok$distanceToTSS <= 0 &
             fok$distanceToTSS > -2000] <- "promoter"
return(fok)
}
