library(poliscidata)
nes <- nes

nes_mod <- glm(voted2012 ~ dem + dem_age6 + dem_edugroup + dem_raceeth,
               family = binomial(link = logit), data = nes)
summary(nes_mod)

# Aggregating for Statistical Tests
nes$voted2012_num <- as.numeric(nes$voted2012) - 1

nes_agg <- aggregate(voted2012_num ~ dem + dem_age6 + dem_edugroup + dem_raceeth,
                     data = nes,
                     FUN = function(x) c(success = sum(x), n = length(x)))

nes_agg$success <- nes_agg$voted2012_num[, "success"]
nes_agg$n <- nes_agg$voted2012_num[, "n"]
nes_agg$failure <- nes_agg$n - nes_agg$success

nes_agg_mod <- glm(cbind(success, failure) ~ dem + dem_age6 + dem_edugroup + dem_raceeth,
                   family = binomial(link = logit),
                   data = nes_agg)

## Residual Deviance Test
# Poor fit if p < 0.05
pchisq(deviance(nes_agg_mod), df.residual(nes_agg_mod), lower.tail = F)

## Pearson Statistic
# Poor fit if p < 0.05
pearson_chi2 <- sum(residuals(nes_agg_mod, type = "pearson")^2)
pearson_chi2

pchisq(pearson_chi2, df.residual(nes_agg_mod), lower.tail = F)

## Hosmer-Lemeshow Test
hoslem.test(nes_mod$y, fitted(nes_mod), g = 10)

## Osius-Rojek Test
osius_rojek(nes_mod)

## Stukel's Test
stukels_test(nes_mod)
