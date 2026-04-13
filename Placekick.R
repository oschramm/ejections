# Placekick
placekick <- read.csv("Placekick.csv")
head(placekick)
tail(placekick)


# Fitting the model
mod_fit <- glm(formula = good ~ distance, family  = binomial(link = logit),
               data = placekick)

summary(mod_fit)  # summary(mod.fit) works too

# Distance and change
mod_fit2 <- glm(formula = good ~ change + distance, family = binomial(link = logit),
                data = placekick)
summary(mod_fit2)
anova(mod_fit2)  # this is sequential (type 1) and so not valid
                     # redone around line 228 (below)
  
# Log-likelihood function
logL <- function(beta, x, Y) {
      pi <- exp(beta[1] + beta[2]*x)/(1+exp(beta[1] + beta[2]*x))
      # Alternatively, could use exp(X%*%beta)/(1+exp(X%*%beta))
      # where X is the design matrix
      sum(Y*log(pi) + (1-Y)*log(1-pi))
    }

############################################################
# Examine the binomial form of the data and re-fit the model

# Find the observed proportion of successes at each distance
w   <- aggregate(x = good ~ distance, data = placekick, FUN = sum)
n   <- aggregate(x = good ~ distance, data = placekick, FUN = length)
w_n <- data.frame(distance = w$distance,
                  success = w$good,
                  trials = n$good,
                  proportion = round(w$good/n$good,4))

head(w_n)
tail(w_n)
sum(w_n$trials)

mod_fit_bin <- glm(formula = proportion ~ distance, weights = trials,
                  family = binomial(link = logit), data = w_n)

pi_hat <- predict(mod_fit_bin, type = "response")
s_res <- rstandard(mod_fit_bin, type = "pearson")
w_n <- data.frame(w_n, pi_hat, s_res)
round(head(w_n), 3)

# Computing the exact lower and upper cumulative prob. of each observed response.

# P (W <= w_m)
prob_less <- pbinom(q = w_n$success, size = w_n$trials, prob = w_n$pi_hat,
                    lower.tail = T)
#P (W >= w_m)
prob_more <- pbinom(q = w_n$success, size = w_n$trials, prob = w_n$pi_hat, 
                    lower.tail = F)

# Minimum of P (W <= w_m) and P (W >= w_m)
w_n$tail_prob <- apply(X = cbind(prob_less, prob_more), MARGIN = 1,
                       FUN = min)
round(head(w_n), 3)

round(w_n[w_n$tail_prob < 0.025,], 3)
round(w_n[abs(w_n$s_res) > 1.96,], 3)

chars <- ifelse(test = w_n$tail_prob < 0.025, yes = 4, no = 1)

# Standardized Pearson residual vs. X plot
plot(x = w_n$distance, y = w_n$s_res, 
     xlab = "Distance",
     ylab = "Standardized Pearson Residuals",
     main = "Standardized Pearson Residuals vs. \n X", pch = chars)
abline(h = c(3,2,0,-2,-3), lty = "dotted", col = "blue")
smooth_stand <- loess(formula = s_res ~ distance, data = w_n, weights = trials)
order_dist <- order(w_n$distance)

lines(x = w_n$distance[order_dist], y = predict(smooth_stand)[order_dist],
      lty = "solid",
      col = "red",
      lwd = 1)

# Standardized Pearson residual vs. pi hat plot
plot(x = w_n$pi_hat, y = w_n$s_res, 
     xlab = "pi",
     ylab = "Standardized Pearson Residuals",
     main = "Standardized Pearson Residuals vs. \n pi", pch = chars)
abline(h = c(3,2,0,-2,-3), lty = "dotted", col = "blue")
smooth_stand <- loess(formula = s_res ~ pi_hat, data = w_n, weights = trials)
order_pi <- order(w_n$pi_hat)

lines(x = w_n$pi_hat[order_dist], y = predict(smooth_stand)[order_dist],
      lty = "solid",
      col = "red",
      lwd = 1)

# Standardized Residuals vs. Linear Predictor
pred <- mod_fit_bin$linear.predictors
w_n$pred <- pred
plot(x = w_n$pred, y = w_n$s_res, 
     xlab = "Linear Predictor",
     ylab = "Standardized Pearson Residuals",
     main = "Standardized Pearson Residuals vs. \n Linear Predictor", pch = chars)
abline(h = c(3,2,0,-2,-3), lty = "dotted", col = "blue")
smooth_stand <- loess(formula = s_res ~ pred, data = w_n, weights = trials)
lines(x = w_n$pred[order_dist], y = predict(smooth_stand)[order_dist],
      lty = "solid",
      col = "red",
      lwd = 1)

mod_fit_bin <- glm(formula = proportion ~ distance, weights = trials,
                   family = binomial(link = logit), data = w_n)

# GOF Statistics
# Deviance/Residual D.F. ratio
rdev <- mod_fit_bin$deviance
dfr <- mod_fit_bin$df.residual
ddf <- rdev/dfr
thresh2 <- 1 + 2*sqrt(2/dfr)
thresh3 <- 1 + 3*sqrt(2/dfr)
c(rdev, dfr, ddf, thresh2, thresh3)
# Doesn't appear to be any serious issues with the ratio

# For GOF tests, we'll have to use packages
# HL Test
library(ResourceSelection)
hoslem.test(w_n$proportion, fitted(mod_fit_bin), g = 10)

# Osius Rojek test and Stukel's Test
library(CLRtools)
# For the Osius-Rojek test, you only need to input your fitted model
osius_rojek(mod_fit_bin)
# Stukel's Test. Only model as input
stukels_test(mod_fit_bin)
