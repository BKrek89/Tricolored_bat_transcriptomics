library(tidyverse)
library(dplyr)
library(emmeans)
library(tidyr)
library(ggplot2)
#sequencing depth by genes discovered
#raw reads by gene number
ggplot(
  dup_graphs_simple,
  aes(x = raw_reads, y = Genes)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 16) +
  labs(
    x = "Raw read count",
    y = "Gene number"
  )

readgenemod <- lm(Genes~raw_reads, dup_graphs_simple)
summary(readgenemod)



#gene number by raw reads and sample group
qc <- expression_data %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  group_by(Sample) %>%
  summarise(
    Sequencing_Depth = sum(Count),
    Genes_1plus = sum(Count > 0),
    Genes_5plus = sum(Count >= 5),
    Genes_10plus = sum(Count >= 10),
    Genes_20plus = sum(Count >= 20)
  )

qc


ggplot(qc, aes(x = Sequencing_Depth, y = Genes_10plus, color = expression_meta$wspec)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "Sequencing depth (mapped reads)",
    y = "Genes detected (≥10 reads)"
  )


ggplot(qc, aes(Sequencing_Depth, Genes_10plus, color = expression_meta$wspec)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18)

summary(lm(Genes_10plus ~ Sequencing_Depth * expression_meta$wspec, data = qc))
confint(lm(Genes_10plus ~ Sequencing_Depth * expression_meta$wspec, data = qc))



####################################################
##total duplication rate plots######################
######################################################
#raw read count x total duplication rate
ggplot(
  dup_graphs_simple,
  aes(x = raw_reads, y = duprate)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 16) +
  labs(
    x = "Raw read count",
    y = "Duplication rate"
  )

countmod <- lm(duprate~raw_reads, dup_graphs_simple)
summary(countmod)

#number of discovered genes by duplication rate
ggplot(
  dup_graphs_simple,
  aes(x = Genes, y = duprate)
) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 16) +
  labs(
    x = "Gene number",
    y = "Duplication rate"
  )

genemod <- lm(duprate~Genes, dup_graphs_simple)
summary(genemod)



#################################################
#per-gene mismatch rate by log2fold change
#################################################

#######################log2fold from post-pesu v pre-mylu
#all genes (post filter)

head(gene_mapping_stats)

#remove two repeat samples of low quality
mismatch <- gene_mapping_stats %>%
  filter(Sample != "1319__515.a1.gen2.bamAligned.sorted.bam")
mismatch1 <- mismatch %>%
  filter(Sample != "1319__516.a1.gen2.bamAligned.sorted.bam")


mismatch2 <- mismatch1 %>%
  left_join(expression_meta, by = "Sample")


mismatch_group <- mismatch2 %>%
  group_by(GeneID, wspec) %>%
  summarise(
    Total_Mismatches = sum(Total_Mismatches, na.rm = TRUE),
    Aligned_Bases = sum(Aligned_Bases, na.rm = TRUE),
    mismatch_rate = Total_Mismatches / Aligned_Bases,
    sex = first(sex),
    batch = first(batch),
    .groups = "drop"
  )


mismatch_group1 <- mismatch_group %>%
  semi_join(resse_pevpremy, by = "GeneID")


plot_data <- mismatch_group1 %>%
  left_join(
    resse_pevpremy %>% select(GeneID, log2FoldChange, padj),
    by = "GeneID"
  )



#absolute value
plot_data <- plot_data %>%
  mutate(abs_log2FC = abs(log2FoldChange))

ggplot(plot_data,
       aes(x = abs_log2FC,
           y = mismatch_rate,
           color = wspec)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_classic(base_size = 18)
model_abs <- lm(
  mismatch_rate ~ abs_log2FC * wspec,
  data = plot_data
)
summary(model_abs)
confint(model_abs)
emtrends(
  model_abs,
  ~ wspec,
  var = "abs_log2FC"
)

pairs(
  emtrends(
    model_abs,
    ~ wspec,
    var = "abs_log2FC"
  )
)

#boxplot of mismatch rate and whether the gene was differentially expressed
plot_data <- plot_data %>%
  mutate(
    DE = abs(log2FoldChange) > 2 & padj < 0.05
  )



plot_data %>%
  group_by(wspec, DE) %>%
  summarise(
    n = n(),
    mean_mismatch = mean(mismatch_rate, na.rm = TRUE),
    median_mismatch = median(mismatch_rate, na.rm = TRUE),
    .groups = "drop"
  )

plot_data <- plot_data %>%
  mutate(
    DE = if_else(
      !is.na(log2FoldChange) & !is.na(padj),
      abs(log2FoldChange) > 2 & padj < 0.05,
      FALSE
    )
  )

ggplot(
  plot_data,
  aes(x = DE, y = mismatch_rate, fill = wspec)
) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  labs(
    x = "Differentially expressed",
    y = "Mismatch rate",
    fill = "Group"
  )


##log2fold from post-pe vs post-my
mismatch2 <- mismatch1 %>%
  left_join(expression_meta, by = "Sample")


mismatch_group <- mismatch2 %>%
  group_by(GeneID, wspec) %>%
  summarise(
    Total_Mismatches = sum(Total_Mismatches, na.rm = TRUE),
    Aligned_Bases = sum(Aligned_Bases, na.rm = TRUE),
    mismatch_rate = Total_Mismatches / Aligned_Bases,
    .groups = "drop"
  )


mismatch_group1 <- mismatch_group %>%
  semi_join(resse_pevpostmy, by = "GeneID")


plot_data <- mismatch_group1 %>%
  left_join(
    resse_pevpostmy %>% select(GeneID, log2FoldChange, padj),
    by = "GeneID"
  )



#absolute value
plot_data <- plot_data %>%
  mutate(abs_log2FC = abs(log2FoldChange))

ggplot(plot_data,
       aes(x = abs_log2FC,
           y = mismatch_rate,
           color = wspec)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_classic()

plot_data <- plot_data %>%
  mutate(
    DE = abs(log2FoldChange) > 2 & padj < 0.05
  )


#boxplot
plot_data <- plot_data %>%
  mutate(
    DE = abs(log2FoldChange) > 2 & padj < 0.05
  )

plot_data %>%
  group_by(wspec, DE) %>%
  summarise(
    n = n(),
    mean_mismatch = mean(mismatch_rate, na.rm = TRUE),
    median_mismatch = median(mismatch_rate, na.rm = TRUE),
    .groups = "drop"
  )

plot_data <- plot_data %>%
  mutate(
    DE = if_else(
      !is.na(log2FoldChange) & !is.na(padj),
      abs(log2FoldChange) > 2 & padj < 0.05,
      FALSE
    )
  )

ggplot(
  plot_data,
  aes(x = DE, y = mismatch_rate, fill = wspec)
) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = "Differentially expressed",
    y = "Mismatch rate",
    fill = "Group"
  )




##log2fold from post-my vs pre-my
mismatch2 <- mismatch1 %>%
  left_join(expression_meta, by = "Sample")


mismatch_group <- mismatch2 %>%
  group_by(GeneID, wspec) %>%
  summarise(
    Total_Mismatches = sum(Total_Mismatches, na.rm = TRUE),
    Aligned_Bases = sum(Aligned_Bases, na.rm = TRUE),
    mismatch_rate = Total_Mismatches / Aligned_Bases,
    .groups = "drop"
  )


mismatch_group1 <- mismatch_group %>%
  semi_join(resse_pomyvpremy, by = "GeneID")


plot_data <- mismatch_group1 %>%
  left_join(
    resse_pomyvpremy %>% select(GeneID, log2FoldChange, padj),
    by = "GeneID"
  )



#absolute value
plot_data <- plot_data %>%
  mutate(abs_log2FC = abs(log2FoldChange))

ggplot(plot_data,
       aes(x = abs_log2FC,
           y = mismatch_rate,
           color = wspec)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_classic()

plot_data <- plot_data %>%
  mutate(
    DE = abs(log2FoldChange) > 2 & padj < 0.05
  )


#boxplot
plot_data <- plot_data %>%
  mutate(
    DE = abs(log2FoldChange) > 2 & padj < 0.05
  )

plot_data %>%
  group_by(wspec, DE) %>%
  summarise(
    n = n(),
    mean_mismatch = mean(mismatch_rate, na.rm = TRUE),
    median_mismatch = median(mismatch_rate, na.rm = TRUE),
    .groups = "drop"
  )

plot_data <- plot_data %>%
  mutate(
    DE = if_else(
      !is.na(log2FoldChange) & !is.na(padj),
      abs(log2FoldChange) > 2 & padj < 0.05,
      FALSE
    )
  )

ggplot(
  plot_data,
  aes(x = DE, y = mismatch_rate, fill = wspec)
) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = "Differentially expressed",
    y = "Mismatch rate",
    fill = "Group"
  )


