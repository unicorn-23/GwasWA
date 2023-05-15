args=commandArgs(T)

data <- read.table(args[1],header = T)
rownum=nrow(data)



limit=args[3]
if(is.na(args[3])){
    limit=0.05/rownum
}

# print(limit)
bon_snp <- subset(data,as.numeric(data$p_wald)<as.numeric(limit))
# bon_snp <- subset(data,as.numeric(data$p_wald)<1e-5)
# bon_snp <- subset(data,as.numeric(data$p_wald)<0.05/rownum)

bon_snp <- bon_snp[order(as.numeric(bon_snp$p_wald)),]
# bon_snp

# e8snp <- subset(e6snp,as.numeric(e6snp$p_wald)<1e-8)
# e8snp <- e8snp[order(as.numeric(e8snp$p_wald)),]

# e8snp

write.table(bon_snp,args[2],row.names=F)
# write.table(e8snp,args[3],row.names=F)
