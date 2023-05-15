args=commandArgs(T)

library("ggplot2")
hwe<-read.table (file=args[1], header=TRUE)


p <- ggplot(data = hwe, mapping = aes(x = -log10(hwe[, 9])))
filter <- -log10(as.numeric(args[2]))

if (max(-log10(hwe[, 9]),na.rm = TRUE) < filter) {
  p + geom_histogram(fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency", x = "HWE P-value(-log10)") 
} else {
  p + geom_histogram( fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency", x = "HWE P-value(-log10)") +
    geom_vline(xintercept = filter, col = "#fa7f6f", linetype = "dashed")
}

ggsave(filename = args[3])
