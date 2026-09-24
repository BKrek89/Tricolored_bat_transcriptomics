
library(BiocManager)
library(dupRadar)
library(ggplot2)

#extract slope and intercepts from dupradar examples
data(dupRadar_examples)

ls()

dm <- dupRadar_examples$dm
dm.bad <- dupRadar_examples$dm.bad

fit_good <- duprateExpFit(dm)
fit_bad  <- duprateExpFit(dm.bad)

example_summary <- data.frame(
  library = c("dupRadar good example",
              "dupRadar bad example"),
  intercept = c(
    fit_good$intercept,
    fit_bad$intercept
  ),
  slope = c(
    fit_good$slope,
    fit_bad$slope
  )
)

example_summary


#graph slope and intercept of RPK vs duplication rate graph
dupRadar_summary$wspec <- factor(dupRadar_summary$wspec,
                                   levels = c("Pre-MYLU", "Post-MYLU", "Post-PESU"))

ggplot(
  dupRadar_summary,
  aes(
    x = intercept,
    y = slope,
    group = wspec
  )
) +
  geom_point(
    aes(color = wspec),
    size = 3
  ) +
  geom_hline(yintercept = 3.186793, linetype = "dashed", colour = "blue") +
  geom_hline(yintercept = 1.442982, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0.04075061, linetype = "dashed", colour = "blue") +
  geom_vline(xintercept = 0.49042274, linetype = "dashed", colour = "red") +
  theme_classic() +
  labs(
    x = "dupRadar intercept",
    y = "dupRadar slope"
  )

#just slope
ggplot(
  dupRadar_summary,
  aes(
    x = wspec,
    y = slope,
    group = wspec
  )
) +
  geom_point(
    aes(color = wspec),
    size = 3
  ) +
  geom_hline(yintercept = 3.186793, linetype = "dashed", colour = "blue") +
  geom_hline(yintercept = 1.442982, linetype = "dashed", colour = "red") +
  theme_classic() +
  labs(
    x = "Sample Group",
    y = "dupRadar slope"
  )


#just intercept
ggplot(
  dupRadar_summary,
  aes(
    x = wspec,
    y = intercept,
    group = wspec
  )
) +
  geom_point(
    aes(color = wspec),
    size = 3
  ) +
  geom_hline(yintercept = 0.04075061, linetype = "dashed", colour = "blue") +
  geom_hline(yintercept = 0.49042274, linetype = "dashed", colour = "red") +
  theme_classic() +
  labs(
    x = "Sample Group",
    y = "dupRadar intercept"
  )

##graph intercept and slope of log(rawcounts) vs duplication rate
dupRadar_rawreg_summary$wspec <- factor(dupRadar_rawreg_summary$wspec,
                                 levels = c("Pre-MYLU", "Post-MYLU", "Post-PESU"))

ggplot(
  dupRadar_rawreg_summary,
  aes(
    x = intercept,
    y = slope,
    group = wspec
  )
) +
  geom_point(
    aes(color = wspec),
    size = 3
  ) +
  theme_classic() +
  labs(
    x = "dupRadar intercept",
    y = "dupRadar slope"
  )
