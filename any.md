'''r
mat <- read.csv("rna_tpm_log2.csv", row.names = 1)
'''

############Generating data frame for visualization#######

'''r
test <- mat["PDK1",]
test <- stack(as.data.frame(test))
test$ind <- colnames(mat)
test$type <- ifelse(test$ind %in% c("MKN1", "Hs746T", "SNU484", "SNU668", "YCC16", "YCC11", "SK4", "SNU1750"), "SEM", "nSEM")
colnames(test) <- c("Expression_Level", "Cell", "Type")
write.csv(test, "PDK1.csv")


library(ggplot2)
library(ggpubr)
'''

###################Compare between groups (Type) + wilcox-rank sum test##############
'''r
p <- ggboxplot(test, x = "Type", y = "Expression_Level", ylab = "Expression_Level",fill = "Type", 
               title = "HMGCR") + 
  theme(axis.text.y=element_text(size=30),
        axis.text.x = element_text(angle = , vjust = 0.5, hjust=1, size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, face = "bold"),
        plot.title = element_text(size = 40, face = "bold"))
a <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               plot.background = element_rect(fill = "transparent", color = "white"),
               panel.background = element_rect(fill = "transparent", color = "white"))

a + stat_compare_means( aes(label = ..p.signif..),) + stat_compare_means(size = 5)
'''

'''r
##################Visualization of all elements - In this case, Descending Order######################
p <- ggplot(test) + 
  geom_boxplot(aes(x= reorder(Cell, -Expression_Level), y = Expression_Level), ) +
  ggtitle("HMGCR") +
  theme(axis.text.y=element_text(size=30),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, face = "bold"),
        plot.title = element_text(size = 40, face = "bold"))

a <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               plot.background = element_rect(fill = "transparent", color = "white"),
               panel.background = element_rect(fill = "transparent", color = "white"))
'''
               
