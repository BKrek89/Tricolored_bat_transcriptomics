library(ggplot2)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(writexl)
library(eulerr)



#create DESeq object
exp_data <- exp_c30_rbad_filt_out

rownames(exp_data) <- exp_data$GeneID
expdatase <- exp_data[, -1]
expdatase <- as.matrix(expdatase)

ddsse <- DESeqDataSetFromMatrix(countData = expdatase, 
                                colData = meta_table_ba, 
                                design = ~ wspec + batch + sex)
#run DESeq analysis
ddsse <- DESeq(ddsse)
resultsNames(ddsse)
resse <- results(ddsse, name="wspec_Post.PESU_vs_Post.MYLU")
resse <- as.data.frame(resse)
write_xlsx(resse, "c:/path/DE/resse.dr.pevpostmy.xlsx")

#rerun above with pre-MYLU as reference
meta_table_ba$wspec <- factor(meta_table_ba$wspec,
                      levels = c("Pre-MYLU", "Post-MYLU", "Post-PESU"))

ddsse <- DESeqDataSetFromMatrix(countData = expdatase, 
                                colData = meta_table_ba, 
                                design = ~ wspec + batch + sex)
ddsse <- DESeq(ddsse)
resultsNames(ddsse)
resse1 <- results(ddsse, name="wspec_Post.PESU_vs_Pre.MYLU")
resse2 <- results(ddsse, name="wspec_Post.MYLU_vs_Pre.MYLU")
resse1 <- as.data.frame(resse1)
resse2 <- as.data.frame(resse2)
write_xlsx(resse2, "path/DE/resse.dr.pomyvpremy.xlsx")
write_xlsx(resse1, "path/DE/resse.dr.pevpremy.xlsx")


#extract significantly differentially expressed genes from each comparison
res_se_filt1 <- resse_dr_pevpostmy[resse_dr_pevpostmy$log2FoldChange > 2 | resse_dr_pevpostmy$log2FoldChange < -2, ]
res_se_filt1.sort <- res_se_filt1[order(res_se_filt1$padj), ]
res.se.sort.filt2 <- filter(res_se_filt1.sort, padj <= 0.05)
nrow(res.se.sort.filt2)
write_xlsx(res.se.sort.filt2, "path/DE/sig.pevpostmy.dr.xlsx")

res_se_filt1 <- resse_dr_pevpremy[resse_dr_pevpremy$log2FoldChange > 2 | resse_dr_pevpremy$log2FoldChange < -2, ]
res_se_filt1.sort <- res_se_filt1[order(res_se_filt1$padj), ]
res.se.sort.filt2 <- filter(res_se_filt1.sort, padj <= 0.05)
nrow(res.se.sort.filt2)
write_xlsx(res.se.sort.filt2, "path/DE/sig.pevpremy.dr.xlsx")

res_se_filt1 <- resse_dr_pomyvpremy[resse_dr_pomyvpremy$log2FoldChange > 2 | resse_dr_pomyvpremy$log2FoldChange < -2, ]
res_se_filt1.sort <- res_se_filt1[order(res_se_filt1$padj), ]
res.se.sort.filt2 <- filter(res_se_filt1.sort, padj <= 0.05)
nrow(res.se.sort.filt2)
write_xlsx(res.se.sort.filt2, "path/DE/sig.pomyvpremy.dr.xlsx")


#######find overlaps of significantly differentially expressed genes among comparisons
#separate each comparison into significantly up and down regulated genes
pevpostmy.up <- filter(sig_pevpostmy_dr, log2FoldChange > 0 )
pevpostmy.do <- filter(sig_pevpostmy_dr, log2FoldChange < 0 )

pevpremy.up <- filter(sig_pevpremy_dr, log2FoldChange > 0)
pevpremy.do <- filter(sig_pevpremy_dr, log2FoldChange < 0)

pomyvpremy.up <- filter(sig_pomyvpremy_dr, log2FoldChange > 0)
pomyvpremy.do <- filter(sig_pomyvpremy_dr, log2FoldChange < 0)

###putative shared response
#extract genes DE in both comparisons to pre-WNS MYLU
join_shr.up <- inner_join(pevpremy.up, pomyvpremy.up, by = "GeneID") %>%
  select(
    GeneID
  )

join_shr.do <- inner_join(pevpremy.do, pomyvpremy.do, by = "GeneID") %>%
  select(
    GeneID
  )
#remove genes that were also DE in pevpostmy comparison
join_shr.up.no <- anti_join(join_shr.up, pevpostmy.do)
join_shr.do.no <- anti_join(join_shr.do, pevpostmy.up)

write_xlsx(join_shr.up, "path/DE/shres.up.dr.xlsx")
write_xlsx(join_shr.do, "path/DE/shres.down.dr.xlsx")

###consistent species differences
join_spec.up <- inner_join(pevpremy.up, pevpostmy.up, by = "GeneID") %>%
  select(
    GeneID
  )

join_spec.do <- inner_join(pevpremy.do, pevpostmy.do, by = "GeneID") %>%
  select(
    GeneID
  )

join_spec.up.no2 <- anti_join(join_spec.up, pomyvpremy.up)
join_spec.do.no2 <- anti_join(join_spec.do, pomyvpremy.do)

write_xlsx(join_spec.up.no, "path/DE/spec.up.dr.xlsx")
write_xlsx(join_spec.do.no, "path/DE/spec.down.dr.xlsx")

###other group
join_other.up <- inner_join(pomyvpremy.up, pevpostmy.do, by = "GeneID") %>%
  select(
    GeneID
  )

join_other.do <- inner_join(pomyvpremy.do, pevpostmy.up, by = "GeneID") %>%
  select(
    GeneID
  )

join_other.up.no <- anti_join(join_other.up, pevpremy.up)
join_other.do.no <- anti_join(join_other.do, pevpremy.do)

write_xlsx(join_spec.up.no, "path/DE/spec.up.dr.xlsx")
write_xlsx(join_spec.do.no, "path/DE/spec.down.dr.xlsx")



join_all.up <- inner_join(join_shr.up, pevpostmy.do, by = "GeneID") %>%
  select(
    GeneID
  )

join_all.do <- inner_join(join_shr.do, pevpostmy.up, by = "GeneID") %>%
  select(
    GeneID
  )




#add functions
shres_down_dr.func <- inner_join(shres_down_dr, Supplementary_Tables_Fin2, by = "GeneID") %>%
  select(GeneID, 
         geneName,
         `function`
         )

write_xlsx(shres_down_dr.func, "path/DE/shres.down.dr.func.xlsx")


##venn diagram
#put number of genes in each group into an excel file and imported
counts.up <- filter(genelist_overlap_dr, reg == "up")
fit <- euler(
  c(
    "A" = counts.up$A,
    "B" = counts.up$B,
    "C" = counts.up$C,
    "A&B" = counts.up$`A&B`,
    "A&C" = counts.up$`A&C`,
    "B&C" = counts.up$`B&C`,
    "A&B&C" = counts.up$`A&B&C`
  ),
  shape = "circle"
)



plot(
  fit,
  fills = list(
    fill = c(
      "A" = "gray70",
      "B" = "gray70",
      "C" = "gray70",
      "A&B" = "#1F78B4",
      "A&C" = "gray70",
      "B&C" = "#1F78B4",
      "A&B&C" = "gray70"
    ),
    alpha = 0.6,
    mode = "disjoint"
  ),
  labels = list(
    labels = c(
      "Post-MYLU vs Pre-MYLU",
      "Post-PESU vs Pre-MYLU",
      "Post-PESU vs Post-MYLU"
    ),
    font = 2,
    cex = 1.2
  ),
  quantities = list(cex = 1.5)
)



counts.do <- filter(genelist_overlap_dr, reg == "do")


fit <- euler(
  c(
    "A" = counts.do$A,
    "B" = counts.do$B,
    "C" = counts.do$C,
    "A&B" = counts.do$`A&B`,
    "A&C" = counts.do$`A&C`,
    "B&C" = counts.do$`B&C`,
    "A&B&C" = counts.do$`A&B&C`
  ),
  shape = "ellipse"
)



plot(
  fit,
  fills = list(
    fill = c(
      "A" = "gray70",
      "B" = "gray70",
      "C" = "gray70",
      "A&B" = "#1F78B4",
      "A&C" = "gray70",
      "B&C" = "#1F78B4",
      "A&B&C" = "gray70"
    ),
    alpha = 0.6,
    mode = "disjoint"
  ),
  labels = list(
    labels = c(
      "Post-MYLU vs Pre-MYLU",
      "Post-PESU vs Pre-MYLU",
      "Post-PESU vs Post-MYLU"
    ),
    font = 2,
    cex = 1.2
  ),
  quantities = list(cex = 1.5)
)


counts.all <- filter(genelist_overlap_dr, reg == "all")


fit <- euler(
  c(
    "A" = counts.all$A,
    "B" = counts.all$B,
    "C" = counts.all$C,
    "A&B" = counts.all$`A&B`,
    "A&C" = counts.all$`A&C`,
    "B&C" = counts.all$`B&C`,
    "A&B&C" = counts.all$`A&B&C`
  ),
  shape = "ellipse"
)



plot(
  fit,
  fills = list(
    fill = c(
      "A" = "gray70",
      "B" = "gray70",
      "C" = "gray70",
      "A&B" = "#1F78B4",
      "A&C" = "gray70",
      "B&C" = "#1F78B4",
      "A&B&C" = "gray70"
    ),
    alpha = 0.6,
    mode = "disjoint"
  ),
  labels = list(
    labels = c(
      "Post-MYLU vs Pre-MYLU",
      "Post-PESU vs Pre-MYLU",
      "Post-PESU vs Post-MYLU"
    ),
    font = 2,
    cex = 1.2
  ),
  quantities = list(cex = 1.5)
)