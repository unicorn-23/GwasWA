args=commandArgs(T)
library(ggplot2)
snpmiss <- read.table(file = args[1], header = TRUE)

p <- ggplot(data = snpmiss, mapping = aes(x = snpmiss[, 5]))
filter <- as.numeric(args[2])
bwidth <- round(max(snpmiss[, 5]) / 50, 5)

if (max(snpmiss[, 5],na.rm = TRUE) < filter) {
  p + geom_histogram(binwidth = bwidth, fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency(log10)", x = "SNP Missingness") + scale_y_log10()
} else {
  p + geom_histogram(binwidth = bwidth, fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency(log10)", x = "SNP Missingness") +
    geom_vline(xintercept = filter, col = "#fa7f6f", linetype = "dashed") + scale_y_log10()
}
ggsave(filename = args[3])
