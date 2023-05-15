args=commandArgs(T)


# 给特征向量加上表头,_eigenvec.xls
eigenvec <- read.table(args[1], quote='\"', comment.char="")
colnames(eigenvec)<-c("FID","sample",paste0("PC",1:(ncol(eigenvec)-2)))
write.table(eigenvec[2:ncol(eigenvec)],file = args[2],sep = "\t",row.names = F,col.names = T,quote = F)

# head(eigenvec)
#    FID sample        PC1        PC2        PC3        PC4        PC5        PC6         PC7        PC8         PC9        PC10
# 1 5837   5837 -0.0240310  0.0323141 -0.0915746  0.1879350 -0.0719426 -0.0901024 -0.00422162 -0.0732364 -0.00215649  0.10340600
# 2 6009   6009  0.0468005 -0.3189510  0.1804640 -0.0347772 -0.1658620 -0.0396738  0.16226600 -0.0461995 -0.07477590  0.00961131
# 3 6898   6898 -0.0309631  0.0595395  0.0295257 -0.0398643  0.0102778 -0.0196374 -0.00467770  0.0330698  0.01581860 -0.12614000
# 4 6900   6900  0.0285349 -0.2342280  0.1719810  0.0573292  0.0956126  0.0281939 -0.10813900  0.0126364  0.59737100  0.03624740
# 5 6901   6901  0.0285098 -0.2341820  0.1719650  0.0573269  0.0955956  0.0282452 -0.10818800  0.0127015  0.59739800  0.03622930
# 6 6903   6903 -0.0166166  0.0285628 -0.0924636  0.2060710 -0.0648266 -0.0810800  0.00199766 -0.0805339  0.00703989  0.10163400

# 特征值
eigenval <- read.table(args[3], quote="\"", comment.char="")
pcs<-paste0("pc",1:nrow(eigenval))
# 计算解释率
percentage<-eigenval$V1/sum(eigenval$V1)*100
eigenval_df<-as.data.frame(cbind(pcs,eigenval[,1],percentage),stringsAsFactors = F)
names(eigenval_df)<-c("pcs","variance","proportation")
eigenval_df$variance<-as.numeric(eigenval_df$variance)
eigenval_df$proportation<-as.numeric(eigenval_df$proportation)
eigenval_df

write.table(eigenval_df,file = args[4],sep = "\t",quote = F,row.names = F,col.names = T)





