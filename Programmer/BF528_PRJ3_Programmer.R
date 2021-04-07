if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
​
BiocManager::install("DESeq2")
library(dplyr)
library(tidyverse)
library(DESeq2)
BiocManager::install("apeglm")
​
#read in each file
SRR1177981<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1177981.txt", header = TRUE)
SRR1177982<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1177982.txt", header = TRUE)
SRR1177983<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1177983.txt", header = TRUE)
SRR1178008<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1178008.txt", header = TRUE)
SRR1178009<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1178009.txt", header = TRUE)
SRR1178010<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1178010.txt", header = TRUE)
SRR1178014<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1178014.txt", header = TRUE)
SRR1178021<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1178021.txt", header = TRUE)
SRR1178047<-read.table("/projectnb/bf528/users/frazzled/project_3/programmer/fCountsSRR1178047.txt", header = TRUE)
​
#combine each file 
sample_file <- full_join(SRR1177981,SRR1177982, by = "Geneid")  #use full_join to join all the txt files based on Geneid 
sample_file2 <-full_join(sample_file,SRR1177983, by = "Geneid")
sample_file3<- full_join(sample_file2,SRR1178008, by = "Geneid")
sample_file4<- full_join(sample_file3,SRR1178009, by = "Geneid")
sample_file5 <- full_join(sample_file4,SRR1178010, by="Geneid")
sample_file6<- full_join(sample_file5,SRR1178014, by="Geneid")
sample_file7<- full_join(sample_file6,SRR1178021, by="Geneid")
samples<-full_join(sample_file7,SRR1178047,by="Geneid") %>% #keep only these columns
  select(Geneid,X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1177981_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1177982_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1177983_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1178008_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1178009_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1178010_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1178014_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1178021_Aligned.sortedByCoord.out.bam,
         X.projectnb.bf528.users.frazzled.project_3.data_curator.STAR_results.SRR1178047_Aligned.sortedByCoord.out.bam)
#rename the columns 
colnames(samples)<-c("Geneid","SRR1177981","SRR1177982","SRR1177983","SRR1178008","SRR1178009","SRR1178010","SRR1178014","SRR1178021","SRR1178047")
samples<-as_tibble(samples)
#export as csv file
write.csv(samples,"SRR_samples.csv")
​
control_sample <- read.csv("/project/bf528/project_3/samples/control_counts.csv")
​
#Pull out each control/sample and merge 
AhR_controls <-subset(control_sample,select = c("Geneid","SRR1178050","SRR1178061","SRR1178063")) 
AhR<-subset(samples,select = c("Geneid","SRR1178008","SRR1178009","SRR1178010")) %>%
  full_join(AhR_controls,by="Geneid")
AhR<-as_tibble(AhR) %>% #convert into tibble df and remove any 0 from the control/samples
  subset(rowSums(AhR==0)==0)
​
#turn geneid col into rownames
AhR<-remove_rownames(AhR) %>%
  column_to_rownames(var="Geneid") 
​
CAR_PXR_controls<-subset(control_sample,select = c("Geneid","SRR1178050","SRR1178061","SRR1178063"))
CAR_PXR <-subset(samples,select = c("Geneid","SRR1178014","SRR1178021","SRR1178047")) %>%
  full_join(CAR_PXR_controls,by="Geneid")
CAR_PXR<-as_tibble(CAR_PXR) %>% #convert into tibble df and remove any 0 from the control/samples
  subset(rowSums(CAR_PXR==0)==0)
​
#turn geneid col into rownames
CAR_PXR<-remove_rownames(CAR_PXR) %>%
  column_to_rownames(var="Geneid") 
​
DNA_DAMAGE_controls<-subset(control_sample,select = c("Geneid","SRR1178004","SRR1178006","SRR1178013"))
DNA_DAMAGE<-subset(samples,select = c("Geneid","SRR1177981","SRR1177982","SRR1177983")) %>%
  full_join(DNA_DAMAGE_controls,by="Geneid")
DNA_DAMAGE<-as_tibble(DNA_DAMAGE) %>% #convert into tibble df and remove any 0 from the control/samples
  subset(rowSums(DNA_DAMAGE==0)==0)
​
#turn geneid col into rownames
DNA_DAMAGE<-remove_rownames(DNA_DAMAGE) %>%
  column_to_rownames(var="Geneid") 
​
#read in metadata file toxgroup 3 and keep only the 
TOX3_META <- read.csv("/project/bf528/project_3/toxgroups/toxgroup_3_rna_info.csv")
#TOX3_META$mode_of_action<-factor(TOX3_META$mode_of_action)
AhR_info <- TOX3_META[1:3,] %>%
  rbind(TOX3_META[10:12,])
AhR_info$mode_of_action<-factor(AhR_info$mode_of_action)
​
​
CAR_PXR_info<- TOX3_META[4:6,] %>%
  rbind(TOX3_META[10:12,])
CAR_PXR_info$mode_of_action<-factor(CAR_PXR_info$mode_of_action)
​
DNA_DAMAGE_info<- TOX3_META[7:9,] %>%
  rbind(TOX3_META[13:15,])
DNA_DAMAGE_info$mode_of_action<-factor(DNA_DAMAGE_info$mode_of_action)
​
#create DESeq object for AhR 
AhR_dds <- DESeqDataSetFromMatrix(
  countData = AhR,
  colData = AhR_info,
  design= ~ mode_of_action
)
# relevel mode_of_action as factor
AhR_dds$mode_of_action <- relevel(AhR_dds$mode_of_action, ref='Control')
# run DESeq
AhR_dds <- DESeq(AhR_dds)
AhR_res <- results(AhR_dds, contrast=c('mode_of_action','AhR','Control'))
AhR_res <- lfcShrink(AhR_dds, coef=2)
AhR_result <- as_tibble(AhR_res,rownames=NA)
AhR_result<- AhR_result %>%
  arrange(padj)
#write csv AhR DESeq results 
write.csv(AhR_result,'AhR_DESeq_results.csv')
#adjusted p values < 0.05
AhR_DESeq <- read.csv("AhR_DESeq_results.csv")
AhR_padj<-AhR_DESeq[AhR_DESeq$padj < 0.05,]
write.csv(AhR_padj, "AhR_DESeq_results_padj.csv")
​
#write csv AhR DESeq normalized counts 
write.csv(counts(AhR_dds,normalized=TRUE),'AhR_norm_DESeq_counts.csv')
​
#create DESeq object for CAR/PXR 
CAR_PXR_dds <- DESeqDataSetFromMatrix(
  countData = CAR_PXR,
  colData = CAR_PXR_info,
  design= ~ mode_of_action
)
​
# relevel mode_of_action as factor
CAR_PXR_dds$mode_of_action <- relevel(CAR_PXR_dds$mode_of_action, ref='Control')
​
#run DESeq 
CAR_PXR_dds <- DESeq(CAR_PXR_dds)
CAR_PXR_res <- results(CAR_PXR_dds, contrast=c('mode_of_action','CAR/PXR','Control'))
CAR_PXR_res <- lfcShrink(CAR_PXR_dds, coef=2)
CAR_PXR_result <- as_tibble(CAR_PXR_res,rownames=NA)
CAR_PXR_result<- CAR_PXR_result %>%
  arrange(padj)
​
#write csv CAR_PXR DESeq results 
write.csv(CAR_PXR_result,'CAR_PXR_DESeq_results.csv')
​
#adjusted p values < 0.05
CAR_PXR_DESeq <- read.csv("CAR_PXR_DESeq_results.csv")
CAR_PXR_padj<-CAR_PXR_DESeq[CAR_PXR_DESeq$padj < 0.05,]
write.csv(CAR_PXR_padj, "CAR_PXR_DESeq_results_padj.csv")
​
#write csv CAR_PXR DESeq normalized counts 
write.csv(counts(CAR_PXR_dds,normalized=TRUE),'CAR_PXR_norm_DESeq_counts.csv')
​
​
#create DESeq object for DNA_DAMAGE
DNA_DAMAGE_dds <- DESeqDataSetFromMatrix(
  countData = DNA_DAMAGE,
  colData = DNA_DAMAGE_info,
  design= ~ mode_of_action
)
​
# relevel mode_of_action as factor
DNA_DAMAGE_dds$mode_of_action <- relevel(DNA_DAMAGE_dds$mode_of_action, ref='Control')
​
#run DESeq 
DNA_DAMAGE_dds <- DESeq(DNA_DAMAGE_dds)
DNA_DAMAGE_res <- results(DNA_DAMAGE_dds, contrast=c('mode_of_action','DNA_Damage','Control'))
DNA_DAMAGE_res <- lfcShrink(DNA_DAMAGE_dds, coef=2)
DNA_DAMAGE_result <- as_tibble(DNA_DAMAGE_res,rownames=NA)
DNA_DAMAGE_result<- DNA_DAMAGE_result %>%
  arrange(padj)
​
#write csv DNA_DAMAGE DESeq results 
write.csv(DNA_DAMAGE_result,'DNA_DAMAGE_DESeq_results.csv')
​
#adjusted p values < 0.05
DNA_DAMAGE_DESeq <- read.csv("DNA_DAMAGE_DESeq_results.csv")
DNA_DAMAGE_padj<-DNA_DAMAGE_DESeq[DNA_DAMAGE_DESeq$padj < 0.05,]
write.csv(DNA_DAMAGE_padj, "DNA_DAMAGE_DESeq_results_padj.csv")
​
#write csv DNA_DAMAGE DESeq normalized counts 
write.csv(counts(DNA_DAMAGE_dds,normalized=TRUE),'DNA_DAMAGE_norm_DESeq_counts.csv')
#DNA_D_ncounts <- read.csv("DNA_DAMAGE_norm_DESeq_counts.csv")
​
#subset top 10 genes arranged by pvalues 
AhR_T10 <- AhR_padj[1:10,]
colnames(AhR_T10)<-c("Geneid","baseMean","log2FoldChange","lfcSE","pvalue","padj")
AhR_T10 <- remove_rownames(AhR_T10) %>%
  column_to_rownames(var="Geneid")
​
CAR_PXR_T10 <- CAR_PXR_padj[1:10,]
colnames(CAR_PXR_T10)<-c("Geneid","baseMean","log2FoldChange","lfcSE","pvalue","padj")
CAR_PXR_T10 <- remove_rownames(CAR_PXR_T10) %>%
  column_to_rownames(var="Geneid")
​
DNA_DAMAGE_T10 <- DNA_DAMAGE_padj[1:10,]
colnames(DNA_DAMAGE_T10)<-c("Geneid","baseMean","log2FoldChange","lfcSE","pvalue","padj")
DNA_DAMAGE_T10 <- remove_rownames(DNA_DAMAGE_T10) %>%
  column_to_rownames(var="Geneid")
