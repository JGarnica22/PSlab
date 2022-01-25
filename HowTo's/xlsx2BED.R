#### load xlsx files with granges dataframe format and convert it to .bed
# set the path where are your xlsx to convert to bed
path <- "~/Bioinformatics/SANTAMARIA_28/out"
#change working directory
setwd(path)
#create directory for bed files
dir.create("bed_files")
# Load xlsx
list.xl <- list.files(path, pattern = ".xlsx")

for (i in list.xl){
  x <- readxl::read_xlsx(i)
  na <- grep("^seqnames|^chr|^start|^end|^stop", names(x), ignore.case = T)
  # or just select the first 3 cols
  # x <- x[,1:3]
  if (length(na)== 3){
  x <- x[,na]
  # generate bed file
  write.table(x, paste0("bed_files/",
                        strsplit(i, split = ".", fixed=T)[[1]][1],
                        ".bed"), row.names = F, col.names = F,
                        quote = F,
                        )
  } else {
    print(paste0("Warning! Error converting ", i ,
                 ". Check the column names in the excel file."))
  }
}
