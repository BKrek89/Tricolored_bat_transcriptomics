library(dplyr)
library(ggplot2)
library(ggbiplot)
library(tidyr)
library(writexl)





#remove low expression. filter row sums
keep <- rowSums(raw_table[,2:56]) >= 55
exp.c99 <- raw_table[keep,]
#remove low expression. filter reads that do not have counts > 10 in at least 5 samples
exp.c55.filt <- exp.c99[apply(exp.c99[, 2:56], 1, function(x) sum(x >= 10) >= 5), ]


library(writexl)
write_xlsx(exp.c55.filt, "path/exp.c55.filt.out.xlsx")


#PCA. add back gene ids. Make sure the first row isn't column names
#normalize with DEseq
library(DESeq2)
rownames(exp_c55_filt_lrs) <- exp_c55_filt_lrs$GeneID
expdata <- exp_c55_filt_lrs[, -1]
expdata <- as.matrix(expdata)
coldata <- meta_table
all(rownames(coldata) == colnames(expdata))

dds <- DESeqDataSetFromMatrix(countData=expdata, 
                              colData=coldata, 
                              design=~wns)
dds <- estimateSizeFactors(dds)
exp.c55.filt.norm <- counts(dds, normalized=TRUE)
exp.c55.filt.norm <- as.data.frame(exp.c55.filt.norm)
library(writexl)
write_xlsx(exp.c55.filt.norm, "path/exp.c55.filt.lrs.norm.xlsx")

#PCA. transposed in excel before being brought back in. 
meta_table$wspec <- factor(meta_table$wspec, levels = c("Pre-MYLU", "Post-MYLU", "Post-PESU"))
exppca <- prcomp(exp.c55.filt.lrs.norm.tran[, 2:5441])



#Biplots

ggbiplot(exppca, 
         ellipse = TRUE, 
         ellipse.alpha = 0,
         groups = meta_table$wspec, 
         shape = meta_table$wns, 
         varname.size = 0, 
         var.axes = FALSE) +
  p$layers[[1]]$aes_params$size <- 4 +
  scale_color_manual(values = c("Pre-MYLU" = "lightblue", 
                                "Post-MYLU" = "royalblue", 
                                "Post-PESU" = "green")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),  
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



p <- ggbiplot(exppca, 
         ellipse = TRUE, 
         ellipse.alpha = 0,
         groups = meta_table$wspec, 
         shape = meta_table$wns, 
         varname.size = 0, 
         var.axes = FALSE) 
  p$layers[[1]]$aes_params$size <- 4 
  p+ scale_color_manual(values = c("Pre-MYLU" = "lightblue", 
                                "Post-MYLU" = "royalblue", 
                                "Post-PESU" = "green")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),  
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



  p <- ggbiplot(
    exppca,
    ellipse = TRUE,
    groups = meta_table$wspec
  )
  
  p$layers


#extract the values of the first 4 principal components
exppcs <- as.data.frame(exppca$x)
exppcs.sc <- scale(exppcs, center = TRUE, scale = TRUE)
exppcs.sc.de <- as.data.frame(exppcs.sc)
exppc4.sc <- exppcs.sc.de[, 1:4]
write_xlsx(exppc4.sc, "path/pcreg/expc4.sc.xlsx")


#pc correlations with variables
corr <- as.data.frame(exppca$rotation)
corr.num <- cbind(" "=rownames(corr), corr)
corr.num5 <- corr.num[, 1:5]
write_xlsx(corr.num5, "path/corr.num.5.xlsx")

# Calculate the variance explained by each principal component
sdev <- exppca$sdev
ve <- sdev^2 / sum(sdev^2)
pc.var <- as.data.frame(ve)
pc.var.num <- pc.var %>%
  mutate(pc = 1:nrow(pc.var))
plot(x=pc.var.num$pc, y=pc.var.num$ve)

write_xlsx(pc.var.num, "c:path/PCreg/pc.var.num.norm.xlsx")








