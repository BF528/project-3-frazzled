rm(list = ls())
#install.packages("BiocManager")
#install.packages("gplots")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("gplots")

#real the normalize csv file 
AhR_norm <- read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/AhR_norm_DESeq_counts.csv")
CAR_PXR_norm<-read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/CAR_PXR_norm_DESeq_counts.csv")
DNA <-read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/DNA_DAMAGE_norm_DESeq_counts.csv")

#real the result csv file 
AHR_result <- read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/AhR_DESeq_results_padj.csv")
CAR_PXR_result<-read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/CAR_PXR_DESeq_results_padj.csv")
DNA_result<-  read.csv("/projectnb/bf528/users/frazzled/project_3/programmer/DNA_DAMAGE_DESeq_results_padj.csv")
row.names(AHR_result)<-AHR_result[,'X']
row.names(CAR_PXR_result)<-CAR_PXR_result[,'X']
row.names(DNA_result)<-DNA_result[,'X']


#merge all the normalize count together 
data_norm <- as.data.frame(full_join(AhR_norm, CAR_PXR_norm, DNA, by="X"), names=1)
data_norm[is.na(data_norm)]<-0
row.names(data_norm)<-data_norm[,1]
data_norm<-data_norm[,-1]

convariantfunc<-function(Cdata){
  N<-nrow(Cdata)
  column_number<-ncol(Cdata)
  logval<-log2(15)
  logpercentofgenes<- (column_number)*0.2
  filter3result<-c();
  counter=0;
  index=1;
  for (x in 1:N){
    std<-sd(Cdata[x,])
    meandata<-mean(as.numeric(Cdata[x,]))
    cvvalue<-std/meandata
    if (cvvalue>0.186){
      filter3result[index]<-x
      index=index+1
    }
    
    
  }
  return(filter3result)
}
finalfilter_3<-convariantfunc(data_norm)
data_norm<-data_norm[finalfilter_3,]



#row.names(data_norm)<-data_norm[,1]
#data_norm<-data_norm[,-1]
plotdata<-data.frame(t(data_norm))
any(is.na(data_norm))
distance<-dist(scale(data.matrix((plotdata))),method ="euclidean")

#cluster 
nest<-hclust(distance, method = "average")
plot(nest)
cluster_cut<-cutree(nest, k=3)
plot(nest)
rect.hclust(nest, k=3)
collection_cutree<-mutate(as.data.frame(plotdata, cluster = cluster_cut))
#Filter 
png('./heatmap.png')


heatmap_matrix<-data.matrix(t(collection_cutree))
colnames(heatmap_matrix)<-rownames(plotdata)

heatmap(heatmap_matrix, main = "Gene expression count heatmap", ylab = "Gene", xlab = "Samples", Rowv = TRUE)

 

dev.off()
#filtering results

#AHR_result, CAR_PXR_result, DNA_result filtering for DAVID enriched pathways

AHR_result <- AHR_result  %>%
filter(padj < 0.05) %>%
arrange(desc(padj))

CAR_PXR_result <- CAR_PXR_result %>%
filter(padj < 0.05) %>%
arrange(desc(padj))

DNA_result <- DNA_result %>%
filter(padj < 0.05) %>%
arrange(desc(padj))

write.csv(row.names(AHR_result), './AhR_newlist.csv')
write.csv(row.names(AHR_result), './CAR_PXR_newlist.csv')
write.csv(row.names(AHR_result), './DNA_newlist.csv')
