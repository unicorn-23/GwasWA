args=commandArgs(T)
library("ggplot2")

tbl=read.table(args[1],header = 1)

hightlight=subset(tbl,Error==min(Error))
p = ggplot(tbl, aes(x = K, y = Error))


p + geom_line() + geom_point(shape=21, color="#82B0D2", fill="#82B0D2", size=3) +
  geom_point(data=hightlight,aes(x=K,y=Error),shape=21, color="#fa7f6f", fill="#fa7f6f", size=4) + labs(y="CV-Error" , x="Groupnum")

ggsave(filename = args[2])
