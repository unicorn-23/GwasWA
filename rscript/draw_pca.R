args=commandArgs(T)
library("ggplot2")
library("colorspace")



eigvec <-read.table(args[1],header = T)
eigval <- read.table(args[2],header = T)
tbl=read.table(args[3])
bestk=args[4]

colnames(tbl) <-c(seq(ncol(tbl) ))
tbl$pop <- apply(tbl, 1, function(t) colnames(tbl)[which.max(t)])

group = paste0("Population",1:bestk)
as.numeric(tbl$pop)
group=group[as.numeric(tbl$pop)]

p=ggplot(data = eigvec, aes(x = PC1, y = PC2, fill = group)) 

p + geom_point(alpha = 1,aes(color=group,shape=group))+guides(fill="none")+
  stat_ellipse(geom = "polygon",alpha=0.5,level = 0.97,aes(fill=group))+
  scale_fill_brewer(palette = "Paired")+scale_color_brewer(palette = "Paired")+
  scale_shape_manual(values=c(15:40))+
  xlab(paste0("PC1(", round(eigval[eigval$pcs == "pc1", 3], 2), "%)")) +
  ylab(paste0("PC2(", round(eigval[eigval$pcs == "pc2", 3], 2), "%)")) +
  # theme_bw()+  theme(panel.grid = element_blank())+
  theme(legend.title = element_blank())

ggsave(filename=args[5])

