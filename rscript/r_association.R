args=commandArgs(T)
library(reshape2)


qfile <- read.table(args[1])

tmp <- read.table(gzfile(args[2]), header = F, stringsAsFactors = F)
ids <- read.table(args[3], header = F, stringsAsFactors = F)

deal_q <- c(ids[1],qfile)
deal_q <- as.data.frame(cbind(ids[1],qfile),stringsAsFactors = F)
header = c("<Trait>",paste0("Q",1:ncol(qfile)))
names(deal_q) <- header
write.table(deal_q,file=args[4],row.names = F, sep = "\t", quote = F)


tmp <- tmp[,c(1,2,4)]
result_matrix <- acast(tmp, V1~V2, value.var="V4", drop = F)
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
result_full <- makeSymm(result_matrix)
diag(result_full) <- 2
result_df <- as.data.frame(result_full)
result_df=cbind(ids[1],result_df)

write.table(result_df, file = args[5], row.names = F, col.names = F, sep = "\t", quote = F)
