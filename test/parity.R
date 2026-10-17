# Optional validation against the downloaded R sources. R is not a runtime dependency.
args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(readxl)
source(file.path(args[1], "R", "classify_asv.R"))
for (folder in list.dirs(args[2], recursive = FALSE, full.names = TRUE)) {
  result <- classify_asv(folder)
  write.table(result, file.path(folder, "r_verdicts.tsv"), sep = "\t",
              row.names = FALSE, quote = TRUE, na = "")
}
source(file.path(args[1], "R", "helpers.R"))
source(file.path(args[1], "R", "hydrophobicity.R"))
write.table(Biostrings::GENETIC_CODE_TABLE[, c("id", "AAs")],
            file.path(args[2], "r_codes.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
data("BLOSUM62", package = "Biostrings")
write.table(BLOSUM62, file.path(args[2], "r_blosum.tsv"), sep = "\t",
            row.names = FALSE, quote = TRUE)
dna <- c("ATGCCAGCA", "TAACCAGCA", "TAAATAAATAA", "ATATGAAGAAGG")
write.table(data.frame(dna = dna, aa = sapply(dna, best_translation, genetic_code = "5")),
            file.path(args[2], "r_translation.tsv"), sep = "\t", row.names = FALSE,
            quote = TRUE, na = "")
aa <- c("ACDEFGHIKLMNPQRSTVWY", "AAAAAAAAA", "LIVMWF")
hydro <- sapply(aa, function(s) get_hydro_metrics(s, hydro_window = 3)$global)
write.table(data.frame(aa = aa, hydro = hydro), file.path(args[2], "r_hydro.tsv"),
            sep = "\t", row.names = FALSE, quote = TRUE)
