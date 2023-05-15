args=commandArgs(T)

#计算杂合度三倍标准差以外的个体
#首先，查看哪些个体在3倍标准差以外，这里是将 所有个体的杂合度平均值 和 
#   群体杂合度的标准差的正负3倍  之和  作为过滤标准， 将群体中个体杂合度离散的个体剔除掉
het <- read.table(args[1], head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-as.numeric(args[3])*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+as.numeric(args[3])*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, args[2], row.names=FALSE)
