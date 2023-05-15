args=commandArgs(T)
library(ggplot2)
library(reshape2)
tbl=read.table(args[1])
tbl2=read.table(args[2])
tbl=cbind(tbl,tbl2)
groupnum=args[3]

colnames(tbl) <-c("strain" ,paste("Population", seq(ncol(tbl)-1 ), sep = ''))
#根据群体含量顺序排列
Population_order = sort(colSums(tbl[,2:ncol(tbl)]), index.return = T, decreasing = T )
tbl <- tbl[,c(1,Population_order$ix + 1 )]

# 第二次排序 根据群体编号再次顺序排列，此时表格的顺序是群体先排在一起，又按照编号排序
tbl$pop <- apply(tbl[,2:ncol(tbl)], 1, function(t) colnames(tbl[,2:ncol(tbl)])[which.max(t)])
tbl <- tbl[with(tbl, order(tbl$pop)),]
tab_group <- tbl

# Determine the factor level according to the current DAT order
tbl$strain <- factor(tbl$strain, levels = tab_group$strain)

data_long <- melt(data = tbl,id.vars=c("strain","pop"),variable.name="Population",value.name="Ancestry")

A=data.frame()
for(i in 1:groupnum){
  A <- rbind(A,c(paste("Population", i, sep = '')))
  
}
colnames(A) <- c("pop")


p <- ggplot(data = data_long, mapping = aes(x = strain  , y = Ancestry, fill = Population))
p + geom_bar(stat = 'identity', position = 'fill',width = 5)+labs(x="",fill = "")+
  
  scale_fill_brewer(palette = "Set3",labels = A$pop) +
  theme(axis.text.x = element_blank())+
  scale_y_continuous(expand = c(0, 0) )
ggsave(filename=args[4],width=20)