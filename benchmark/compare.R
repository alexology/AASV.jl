args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
reference <- args[2]
suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(dplyr)
  library(tibble)
})
for (name in c("helpers.R", "hydrophobicity.R", "calculate_asv_distance.R",
               "asv_functional_structure.R", "classify_asv.R")) {
  source(file.path(reference, "R", name))
}
times <- data.frame(language=character(), operation=character(), n=integer(),
                    phase=character(), seconds=double())
pipeline_only <- length(args)>=3 && args[3]=="--pipeline-only"
if (pipeline_only && file.exists(file.path(root,"r_times.csv"))) {
  previous <- read.csv(file.path(root,"r_times.csv"))
  times <- previous[previous$operation!="functional_pipeline",]
}
measure <- function(operation, n, fn) {
  elapsed <- system.time(fn())[["elapsed"]]
  times <<- rbind(times, data.frame(language="R",operation=operation,n=n,phase="first",seconds=elapsed))
  for (i in seq_len(3)) {
    gc()
    elapsed <- system.time(fn())[["elapsed"]]
    times <<- rbind(times, data.frame(language="R",operation=operation,n=n,phase="warm",seconds=elapsed))
  }
  cat(operation, n, "warm median", median(tail(times$seconds,3)), "\n")
  write.csv(times,file.path(root,"r_times.csv"),row.names=FALSE)
}
if (!pipeline_only) {
  for (n in c(20L,80L)) {
    proteins <- readAAStringSet(file.path(root,paste0("proteins_",n,".fasta")))
    measure("multiple_alignment",n,function() AlignSeqs(proteins,processors=1,verbose=FALSE))
  }
}
asvs <- read.csv(file.path(root,"asvs.csv"),check.names=FALSE)
out <- file.path(root,"r_output")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
measure("functional_pipeline",10L,function() {
  suppressMessages(asv_functional_structure(asvs, output_dir=out,
    alignment_dir=file.path(root,"references"),processors=1))
  classify_asv(out)
})
write.csv(times,file.path(root,"r_times.csv"),row.names=FALSE)
write.csv(classify_asv(out),file.path(root,"r_verdicts.csv"),row.names=FALSE)
