args=commandArgs(T)
library("colorspace")

tab <- read.table(args[1])

group = paste0("G",1:ncol(tab))
color<-qualitative_hcl(as.numeric(ncol(tab)),palette = "Dark 3")
pch<-c(15:(15+ncol(tab)))

# 根据每个样本的最大概率判断属于哪个群体
t <- apply(tab, 1, function(x){which.max(x)})

adinfo<-as.data.frame(cbind(group[t],color[t],pch[t]),stringsAsFactors = F)

names(adinfo)<-c("group","color","pch")
write.csv(adinfo,file=args[2],row.names = F)