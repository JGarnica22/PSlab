# what species use to look for motif into JASPAR database
# generally, most of the motifs are described in humans and are highly
# conserved so, in general, is preferable to use even with mouse data.

species <- "mouse" # say here either "mouse" or "human"
if (species == "mouse") {
  bsgen <- "BSgenome.Mmusculus.UCSC.mm10"
  gm <- "mm10"
} else {
  bsgen <- "BSgenome.Hsapiens.UCSC.hg38"
  gm <- "hg38"
}
bioc.packages <- c(bsgen, "JASPAR2020", "TFBSTools", "motifmatchr")
for (i in bioc.packages) {
  # if (!require(i, character.only = T)) {
  #   BiocManager::install(i)
  #   print(paste(i,"just installed"))
  # } else {
  #   print(paste(i,"was already installed"))
  # }
  library(i, character.only = T)
}

mouse <- 10090
human <- 9606


# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = human , # set species
              all_versions = FALSE))


# let's use motifmatchr to discover what described motif
# are found within a certain peak set.

# This method has two mandatory arguments:
# Position weight matrices or position frequency matrices, stored in the PWMatrix,
# PFMatrix, PWMatrixList, or PFMatrixList objects from the TFBSTools package
# Either a set of genomic ranges (GenomicRanges or RangedSummarizedExperiment object)
# or a set of sequences (either DNAStringSet, DNAString, or simple character vector)
# If the second argument is a set of genomic ranges, a genome sequence is also required.

## load file of peaks to compare and generate granges object
a <- readxl::read_xlsx("out/Diff_Tet_ATAC_Tet_vs_TH0_peaks_annotated.xlsx")
gr <- GRanges(seqnames = a$seqnames, # chromosomes names
              ranges = paste0(a$start,"-", a$end), #start position - end position
              strand = a$strand) # strand, can be NULL


motifs <- matchMotifs(pwm, gr, genome = gm, out = "scores")

# get the name of TF for the matched motifs
for (i in motifs@colData@rownames){
  print(pwm[[i]]@name)
}
