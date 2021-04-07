#Project 3 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggplot2)
library(ggpubr)
#Part 3.5 create boxplots for each sample 
samples_Bplots<- samples %>%
  reshape2::melt(id.vars="Geneid") %>%
  ggplot(aes(x=variable, y=value)) +
  geom_boxplot() +
  scale_y_log10() +
  ggtitle("Toxgroup 3 Read Count Distributions") +
  xlab("Samples") + theme(axis.text.x = element_text(angle = 30,hjust = 1))
samples_Bplots
​
#Part 4.6 create scatter plot
AhR_scatter<- ggplot(AhR_padj,aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point() +
  ggtitle("AhR Differentially Expressed Genes") +
  xlab("Log2 Fold Change") + ylab("-log10 P-adjusted values")
AhR_scatter
​
CAR_PXR_scatter <- ggplot(CAR_PXR_padj,aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point() +
  ggtitle("CAR/PXR Differentially Expressed Genes") + 
  xlab("Log2 Fold Change") + ylab("-log10 P-adjusted values")
CAR_PXR_scatter
​
DNA_DAMAGE_scatter <- ggplot(DNA_DAMAGE_padj, aes(x=log2FoldChange,y=-log10(padj))) +
  geom_point() +
  ggtitle("DNA Damage Differentially Expressed Genes") +
  xlab("Log2 Fold Change") + ylab("-log10 P-adjusted values")
DNA_DAMAGE_scatter
​
#Part 4.6 create histograms for each sample
Ahr_Hist <- ggplot(AhR_padj, aes(AhR_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color = "white", fill= "red") +
  xlab("Log2 Fold Change") + ylab("Values") +
  ggtitle("AhR Fold Change Values")
Ahr_Hist
​
CAR_PXR_Hist <- ggplot(CAR_PXR_padj, aes(CAR_PXR_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color = "white", fill = "blue") +
  xlab("Log2 Fold Change") + ylab("Values") +
  ggtitle("CAR/PXR Fold Change Values")
CAR_PXR_Hist
​
DNA_DAMAGE_Hist <- ggplot(DNA_DAMAGE_padj, aes(DNA_DAMAGE_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color = "white", fill = "dark green") +
  xlab("Log2 Fold Change") + ylab("Values") +
  ggtitle("DNA DAMAGE Fold Change Values")
DNA_DAMAGE_Hist
​
#Combine images into one figure 
Hist <- cowplot::plot_grid(Ahr_Hist,CAR_PXR_Hist,DNA_DAMAGE_Hist, labels = c("A","B","C"))
Hist
​
Scatter <- cowplot::plot_grid(AhR_scatter,CAR_PXR_scatter,DNA_DAMAGE_scatter, labels= c("A","B","C"))
Scatter
