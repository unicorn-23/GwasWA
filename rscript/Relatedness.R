args=commandArgs(T)

png(args[1])
relatedness = read.table(args[2], header=T)
par(pch=16, cex=1)
with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
with(subset(relatedness,RT=="OT") , points(Z0,Z1,col=4))
with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))

legend(1,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=c(4,3))

png(args[3])
relatedness_zoom = read.table(args[4], header=T)
par(pch=16, cex=1)
with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n"))
with(subset(relatedness_zoom,RT=="OT") , points(Z0,Z1,col=4))
with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))
legend(0.02,1, xjust=1, yjust=1, legend=levels(as.factor(relatedness$RT)), pch=16, col=c(4,3))


png(args[5])
relatedness = read.table(args[6], header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")  
dev.off() 

