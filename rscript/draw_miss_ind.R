args = commandArgs(T)
library(ggplot2)
indmiss <- read.table(file = args[1],  header = TRUE)

p <- ggplot(data = indmiss, mapping = aes(x = indmiss[, 6]))
bwidth <- round(max(indmiss[, 6]) / 50, 5)
filter <- as.numeric(args[2])

if (max(indmiss[, 6],na.rm = TRUE) < filter) {
  p + geom_histogram(binwidth = bwidth, fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency(log10)", x = "Individual Missingness") + scale_y_log10()
} else {
  p + geom_histogram(binwidth = bwidth, fill = "#82B0D2", color = "#82B0D2", alpha = 0.8) + labs(y = "Frequency(log10)", x = "Individual Missingness") +
    geom_vline(xintercept = filter, col = "#fa7f6f", linetype = "dashed") + scale_y_log10()
}
ggsave(filename = args[3])
