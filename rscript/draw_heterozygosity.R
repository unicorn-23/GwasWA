args=commandArgs(T)

library("ggplot2")
het <- read.table(args[1], head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."

p <- ggplot(data = het, mapping = aes(x = HET_RATE))

myhet=as.numeric(args[2])

leftfilter <- mean(het$HET_RATE)-myhet*sd(het$HET_RATE)
rightfilter <- mean(het$HET_RATE)+myhet*sd(het$HET_RATE)


if (max(het$HET_RATE ,na.rm = TRUE) < rightfilter && (min(het$HET_RATE ,na.rm = TRUE) > leftfilter)) {
  p + geom_histogram(fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency", x = "Heterozygosity Rate") 
} else {
  p + geom_histogram( fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency", x = "Heterozygosity Rate") +
    geom_vline(xintercept = leftfilter, col = "#fa7f6f", linetype = "dashed")+
    geom_vline(xintercept = rightfilter, col = "#fa7f6f", linetype = "dashed")
}

ggsave(filename =args[3])
