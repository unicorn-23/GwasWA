args=commandArgs(T)

gender <- read.table(args[1], header = T, as.is = T)
library("ggplot2")

p <- ggplot(data = gender, mapping = aes(x = gender[, 6], fill = as.character(PEDSEX),color = as.character(PEDSEX)))
p + geom_histogram(alpha = 0.8)+ labs(fill = "", y = "Count", x = "F-score") + 
  scale_fill_manual(values = c("1" = "#82B0D2", "2" = "#fa7f6f","0"="#999999"), labels = c("Male", "Female","Problem"))+
  scale_color_manual(values = c("1" = "#82B0D2", "2" = "#fa7f6f","0"="#999999"))+guides(color="none")

ggsave(filename = args[2])
