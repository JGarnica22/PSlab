
RPG_list <- list.files(path=paste0(getwd(),"/mouse/alignment_counts"), pattern= "*ReadsPerGene*")

for (g in 1:length(RPG_list)){
  x <- read.table(paste0(getwd(), "/mouse/alignment_counts/", RPG_list[g]),
                  header = T,
                  sep = "\t",
                  quote = "", 
                  dec = ".")
  x1 <- x[-(1:4),c(1,2)]
  names(x1) <- c("Ensembl_id", strsplit(as.character(RPG_list[g]), split="_", fixed=TRUE)[[1]][1])
  if (g == 1){
    all_reads <- x1  
  } else {
    all_reads <- merge(all_reads, x1, by="Ensembl_id")
  }
}

write.table(all_reads, 
            file = paste0(getwd(),"/mouse/alignment_counts/Reads_all_samples.txt"),
            sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
