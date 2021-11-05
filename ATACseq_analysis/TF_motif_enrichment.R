library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)



# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = c(9606,10090) , all_versions = FALSE)
)



# add motif information
mo <- AddMotifs(scdata.f, genome = get(bsgen), pfm = pwm)

for (i in toupper(features[45:98])){
  scdata.f <- Footprint(
    object = scdata.f,
    motif.name = i,
    genome = get(bsgen))
  p2[[u]] <- PlotFootprint(scdata.f,
                           features = i)
  u <- u+1
  
}
