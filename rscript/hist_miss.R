args=commandArgs(T)
library(ggplot2)
indmiss<-read.table(file=args[1], header=TRUE)
snpmiss<-read.table(file=args[2], header=TRUE)


png(args[3]) #indicates png format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness") #selects column 6, names header of file

png(args[4]) 
hist(snpmiss[,5],main="Histogram SNP missingness")  


p <- ggplot(data = indmiss,mapping = aes(x = indmiss[,6]))
if(max(indmiss[,6])<0.02){
  p + geom_histogram(bins=30 , fill="#82B0D2" ,color="#82B0D2",alpha=0.8) + labs(y="Frequency" , x="individual missingness")
}else{
  p + geom_histogram(bins=30 , fill="#82B0D2" ,color="#82B0D2",alpha=0.8) + labs(y="Frequency" , x="individual missingness") + geom_vline(xintercept=0.02,col="#fa7f6f")
}
ggsave(filename=args[3])



p <- ggplot(data = snpmiss,mapping = aes(x = snpmiss[,5]))
if(max(snpmiss[,5])<0.02){
  p + geom_histogram(fill="#82B0D2" ,color="#82B0D2",alpha=0.8) + labs(y="Frequency" , x="SNP missingness") + scale_y_log10()
}else{
  p + geom_histogram(fill="#82B0D2" ,color="#82B0D2",alpha=0.8) + labs(y="Frequency" , x="SNP missingness") + geom_vline(xintercept=0.02,col="#fa7f6f") + scale_y_log10()  
}
ggsave(filename=args[4])
