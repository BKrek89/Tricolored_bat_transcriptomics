library(ggplot2)
library(DHARMa)
library(emmeans)
library(dplyr)



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

