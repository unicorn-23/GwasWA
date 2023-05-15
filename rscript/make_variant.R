args=commandArgs(T)

data <- read.table(args[1],header = T)



variant <- cbind(data$chr,data$ps,data$rs,data$allele0,data$allele1)

if(nrow(variant)>1){
    variant <- data.frame(variant,'.','.','.')
}
write.table(variant,args[2],row.names = F,col.names = F,quote=F)