args=commandArgs(T)


library(ggplot2)
maf_freq <- read.table(args[1], header =TRUE, as.is=T)

p <- ggplot(data = maf_freq, mapping = aes(x = maf_freq[, 5]))

filter <- as.numeric(args[2])
if (max(maf_freq[, 5],na.rm = TRUE) < filter) {
  p + geom_histogram(binwidth = 0.01, fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency(log10)", x = "Minor Allele Frequency") + scale_y_log10()
} else {
  p + geom_histogram(binwidth = 0.01, fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency(log10)", x = "Minor Allele Frequency") + 
    geom_vline(xintercept = filter, col = "#fa7f6f", linetype = "dashed") + scale_y_log10()
}

ggsave(filename = args[3])
