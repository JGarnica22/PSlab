# Extract from a data frame with different sequences a combined FASTA file to use as a intput for ClustalW aligment.

muscle <- c()
m <- 1
TR <- "TRB"
for (i in imgtA[imgtA$TCR == TR,"cell_name"]){
  muscle[m] <- paste0(">", imgtA[imgtA$cell_name == i & 
                                   imgtA$TCR == TR,"cell_name"])
  muscle[m+1] <- imgtA[imgtA$cell_name == i & 
                         imgtA$TCR == TR,"CDR3.IMGT"]
  m <- m+2
}

write.table(muscle, "out/muslce.fasta", row.names = F,
            col.names = F)
            
# Go to ClustalW website and run it!

ma <- read.table("out/clwB.clw",
                    sep = "\t", quote = "",
                    dec = ".", header = T, na.strings = T, fill = T)
ma <- ma %>% separate(CLUSTAL.multiple.sequence.alignment.by.MUSCLE..3.8., sep = "     ",
                into = c("cell_name", "muscle"))
imgtAb <- imgtA %>% filter(TCR == "TRB")
imgtAb <- merge(imgtAb, ma, by = "cell_name", all.y = F)
imgtAb <- imgtAb[!duplicated(imgtAb),]

# Repeat above for the other chain
imgt <- rbind(imgtAa, imgtAb)

writexl::write_xlsx(imgt, "out/imgt_analysis.xlsx", col_names = T)
