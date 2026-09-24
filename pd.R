library(dplyr)
library(ggplot2)
library(ggbiplot)
library(tidyr)
library(writexl)
library(DHARMa)



#remove low expression. filter row sums
keep <- rowSums(pd_raw[,2:54]) >= 53
pd_raw.fil1 <- pd_raw[keep,]

#remove low expression. filter row values
pd_raw.fil2 <- pd_raw.fil1[apply(pd_raw.fil1[, 2:54], 1, function(x) sum(x >= 10) >= 1), ]

#remove any read from post-WNS samples that appeared in pre-WNS samples
pd_fil.nopre <- pd_raw.fil2[!apply(pd_raw.fil2[ , 2:26], 1, function(row) any(row > 0)), ]

write_xlsx(pd_fil.nopre, "pd_fil.nopre.xlsx")




####regression

pdmod <- lm(pdtprop~wspec+Sex+Batch, data_pd)
summary(pdmod)
confint(pdmod)
#Dharma test model fit
simvarg <- simulateResiduals(fittedModel = pdmod, plot = F)
residuals(simvarg)
plot(simvarg)
testDispersion(simvarg)
testOutliers(simvarg)

#normal distribution didn't fit well. Tried Gamma with a log link.
pdgam <- glm(pdtprop~wspec+Sex+Batch, data_pd, family=Gamma(link="log"), na.action = "na.fail")
summary(pdgam)
confint(pdgam)
simvarg <- simulateResiduals(fittedModel = pdgam, plot = F)
residuals(simvarg)
plot(simvarg)
testDispersion(simvarg)
testOutliers(simvarg)

