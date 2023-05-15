args=commandArgs(T)
library("pheatmap")

dat<- read.table(args[1])
png(args[2],width=800,height=800)
pheatmap(dat, fontsize_row = 0.3, fontsize_col = 0.3)
dev.off()
