# Development-only export of the BLOSUM62 data used by the R reference.
# ProfileAlignments does not depend on R at runtime.
data("BLOSUM62", package = "Biostrings")
order <- strsplit("ARNDCQEGHILKMFPSTWYVBZX*", "")[[1]]
destination <- file.path("..", "ProfileAlignments.jl", "data", "BLOSUM62.tsv")
before <- as.matrix(read.table(destination, header = TRUE, row.names = 1, check.names = FALSE))
after <- BLOSUM62[order, order]
changed <- which(before != after, arr.ind = TRUE)
print(data.frame(row = rownames(before)[changed[, 1]],
                 column = colnames(before)[changed[, 2]],
                 old = before[changed], reference = after[changed]))
con <- file(destination, "w")
writeLines(paste(order, collapse = "\t"), con)
for (aa in order) writeLines(paste(c(aa, after[aa, ]), collapse = "\t"), con)
close(con)
