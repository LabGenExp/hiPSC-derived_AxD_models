library(stringr)
library(ggplot2)
library(ggpubr)

### cortical
df <- read.table("cortical_qPCR_log2FC.txt", sep = ",", header = T)

df[c('cond', 'batch')] <- str_split_fixed(df$sample, ' ', 2)
df <- df[,c("NCAN","NNAT","FOXG1","DCX","AQP4","NFIA","OLIG2","TNNT2","cond","batch")]
df$cond <- as.factor(df$cond)
df$cond <- factor(df$cond,c("CTRL","AxD"))

plot_list <- lapply(seq_along(colnames(df)[1:8]), function(x)
  ggplot(df, aes_string("cond", colnames(df)[x])) + 
  geom_boxplot(linewidth = 1) + 
  geom_point(size = 3, aes(color = batch)) +
  ylab(expression("Log"[2]* " normalized relative expression")) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 14),
        axis.title.y = element_blank(), #Log2 normalized relative expression
        axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  ggtitle(colnames(df)[x]) +
  scale_y_continuous(limits = c(NA, 2*max(df[,x]))) +
  stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE), label = "p.signif", label.x = 1.5, size = 6, label.y = 1.1*max(df[,x]), 
                     hide.ns = F, comparisons = list(c("CTRL","AxD")))
)

ggsave("qPCR_cortical.png",
       egg::ggarrange(plots=plot_list, nrow = 2),
       device = "png", width = 8.5, height = 6.5)

### cerebral
df <- read.table("cerebral_qPCR_log2FC.txt", sep = ",", header = T)

df[c('cond', 'batch')] <- str_split_fixed(df$sample, ' ', 2)
df <- df[,c(2:22,1,23,24)]
df <- df[,c("AQP4","TRPM3","PHOX2B","HOPX","MYF5","CLPS","cond","batch")]
df$cond <- as.factor(df$cond)
df$cond <- factor(df$cond,c("CTRL","AxD"))

plot_list <- lapply(seq_along(colnames(df)[1:6]), function(x)
  ggplot(df, aes_string("cond", colnames(df)[x])) + 
    geom_boxplot(linewidth = 1) + 
    geom_point(size = 3, aes(color = batch)) +
    ylab(expression("Log"[2]* " normalized relative expression")) +
    theme_light() +
    theme(axis.title.x = element_blank(),
          #axis.title.y = element_text(size = 14),
          axis.title.y = element_blank(), #Log2 normalized relative expression
          axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1, colour = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    ggtitle(colnames(df)[x]) +
    scale_y_continuous(limits = c(NA, 2*max(df[,x]))) +
    stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE), label = "p.signif", label.x = 1.5, size = 6, label.y = 1.1*max(df[,x]), 
                       hide.ns = F, comparisons = list(c("CTRL","AxD")))
)

ggsave("qPCR_cerebral.png",
       egg::ggarrange(plots=plot_list, nrow = 2),
       device = "png", width = 6.5, height = 6.5)



