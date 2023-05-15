args=commandArgs(T)

library(qqman)  ## 加载qqman包
data <- read.table(args[1], head=TRUE) 
data$snp <- paste("chr",data$chr,"_",data$ps,sep="")

png(args[2], width = 2000,height = 800)
manhattan(data,chr="chr",snp="snp",bp="ps",p="p_wald",col = c("#82B0D2","#FFBE7A"),
          annotatePval = 0.00000015,annotateTop = F,ylim = c(0, 10),
          suggestiveline=F,genomewideline=-log10(1.5e-7)
          )
dev.off()


png(args[3], width = 800,height = 800)
qq(data$p_wald, main = "Q-Q plot", col = "#82B0D2", las = 1)
dev.off()
