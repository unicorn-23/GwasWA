args=commandArgs(T)
library(reshape2)

tmp <- read.table(gzfile(args[1]), header = F, stringsAsFactors = F)
ids <- read.table(args[2], header = F, stringsAsFactors = F)
tmp <- tmp[,c(1,2,4)]
result_matrix <- acast(tmp, V1~V2, value.var="V4", drop = F)
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
result_full <- makeSymm(result_matrix)
diag(result_full) <- 2
result_df <- as.data.frame(result_full)
row.names(result_df) <- ids$V2
colnames(result_df) <- ids$V2
write.table(result_df, file = args[3], row.names = T, col.names = NA, sep = "\t", quote = F)


