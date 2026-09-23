library(ggplot2)
library(DHARMa)
library(emmeans)
library(dplyr)
library(tidyr)


###run models of PC1 and PC2 by sample group with sex and sequencing batch as covariates 
normful <- lm(PC1~wspec+sex+batch, pc_samp, na.action = "na.fail")
summary(normful)
confint(normful)
emmeans(normful, list(pairwise ~ wspec), adjust = "tukey")
confint(emmeans(normful, list(pairwise ~ wspec), adjust = "tukey"))

normful2 <- lm(PC2~wspec+sex+batch, pc_samp, na.action = "na.fail")
summary(normful2)
confint(normful2)
emmeans(normful2, list(pairwise ~ wspec), adjust = "tukey")
confint(emmeans(normful2, list(pairwise ~ wspec), adjust = "tukey"))


#evaluate model fit with DHARMa

simvar1l <- simulateResiduals(fittedModel = normful, plot = F)
residuals(simvar1l)
plot(simvar1l)
testDispersion(simvar1l)
testOutliers(simvar1l)


simvar2 <- simulateResiduals(fittedModel = normful2, plot = F)
residuals(simvar2)
plot(simvar2)
testDispersion(simvar2)
testOutliers(simvar2)



#pc1 model did not fit well with a normal distribution. Tried Gamma with a log link
min(pc_samp$PC1)
pc_samp$PC11 <- pc_samp$PC1+1

gamful <- glm(PC11~wspec+sex+batch, pc_samp, family=Gamma(link="log"), na.action = "na.fail")
summary(gamful)
confint(gamful)
emmeans(gamful, list(pairwise ~ wspec), adjust = "tukey")
confint(emmeans(gamful, list(pairwise ~ wspec), adjust = "tukey"))

#evaluate model fit with DHARMa
simvarg <- simulateResiduals(fittedModel = gamful, plot = F)
residuals(simvarg)
plot(simvarg)
testDispersion(simvarg)
testOutliers(simvarg)




###create plots for the genes most associated with each PC.
#S100 counts per sample extracted manually and imported.
s100_counts$group <- factor(s100_counts$group, levels = c("Pre-MYLU", "Post-MYLU", "Post-PESU"))


# Reshape the data to long format, using gene_name and set factor levels
wns_long <- s100_counts %>%
  pivot_longer(cols = c(S100A8, S100A12, S100A4), 
               names_to = "gene_name", 
               values_to = "value") %>%
  mutate(gene_name = factor(gene_name, levels = c("S100A12", "S100A8", "S100A4")))

# Create the faceted boxplot

ggplot(data = wns_long, aes(x = group, y = value, color = group)) +
  geom_boxplot(linewidth = 1.2) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    legend.position = "none",
    
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 18, face = "bold")
  ) +
  labs(x = NULL, y = "Normalized read count") +
  facet_wrap(~ gene_name, scales = "free_y") +
  scale_color_manual(values = c("lightblue", "royalblue", "green"))
