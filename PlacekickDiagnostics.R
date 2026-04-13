#  Example: Placekicking, diagnostics

options(width = 60)  # Formatting for book - 60 characters per line

placekick <- read.csv("Placekick.csv")
head(placekick)
tail(placekick)

#####################################################################
# Aggregate data according to variables we plan to put into the model

# Find the observed proportion of successes at each distance
w <- aggregate(x = good ~ distance, data = placekick, FUN = sum)
n <- aggregate(x = good ~ distance, data = placekick, FUN = length)
w.n <- data.frame(distance = w$distance, success = w$good, trials = n$good, prop = round(w$good/n$good,4))
head(w.n)
tail(w.n)

# Then refit model using aggregated data
mod.fit.bin <- glm(formula = success/trials ~ distance, weights = trials, family = binomial(link = logit), data = w.n)
summary(mod.fit.bin)

# Show how to find the hat matrix using matrix algebra
X <- model.matrix(mod.fit.bin)
# mod.fit.bin$weights is a shortcut for n*hat(pi)*(1-hat(pi))
# Compare: w.n$trials*mod.fit.bin$fitted.values*(1-mod.fit.bin$fitted.values)
V <- diag(mod.fit.bin$weights)
H <- sqrt(V)%*%X%*%solve(t(X)%*%V%*%X)%*%t(X)%*%sqrt(V)
head(diag(H))
head(hatvalues(mod.fit.bin))  # Matches

#####################################################################
# Create plot with fit overlaid onto confidence intervals for each distance
# Get plotting quantities
pi.hat <- predict(mod.fit.bin, type = "response")
s.res <- rstandard(mod.fit.bin, type = "pearson")
w.n <- data.frame(w.n, pi.hat, s.res)
round(head(w.n), digits = 3)

# Compute exact lower and upper cumulative probability of each observed response
#  (How extreme is the residual?)
# P(W <= w_m)
prob.smaller <- pbinom(q = w.n$success, size = w.n$trials, prob = w.n$pi.hat, lower.tail = TRUE)
# P(W >= w_m) 
prob.higher <- pbinom(q = w.n$success, size = w.n$trials, prob = w.n$pi.hat, lower.tail = FALSE) + dbinom(x = w.n$success, size = w.n$trials, prob = w.n$pi.hat)
# Mininum of P(W <= w_m) and P(W >= w_m) 
w.n$tail.prob <- apply(X = cbind(prob.smaller, prob.higher), MARGIN = 1, FUN = min)


sres.prob <- cbind(pbinom(q = w.n$success, size = w.n$trials, prob = w.n$pi.hat, lower.tail = TRUE),
                  pbinom(q = w.n$success, size = w.n$trials, prob = w.n$pi.hat, lower.tail = FALSE)
                  + dbinom(x = w.n$success, size = w.n$trials, prob = w.n$pi.hat))

w.n <- data.frame(w.n, tail.prob = apply(X = sres.prob, MARGIN = 1, FUN = min))
round(head(w.n), digits = 3)

round(w.n[w.n$tail.prob < 0.025,], digits = 3)
round(w.n[abs(w.n$s.res) > 1.96,], digits = 3)

####################################################################
# Residual plots
# Add linear predictor to w.n for plotting
lin.pred <- mod.fit.bin$linear.predictors
w.n <- data.frame(w.n,lin.pred)
round(head(w.n), digits = 3)


dev.new(height = 7, width = 13, pointsize=20)

# Plotting characters, X for extreme residuals
chars <- ifelse(test = w.n$tail.prob < 0.025, yes = 4, no = 1)
par(mfrow = c(1,3))

# Standardized Pearson residual vs x plot
plot(x = w.n$distance, y = w.n$s.res, xlab = "Distance", ylab = "Standardized Pearson residuals",
   main = "Standardized residuals vs. \n X", pch=chars
   )
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "blue")
smooth.stand <- loess(formula = s.res ~ distance, data = w.n, weights = trials)
# Make sure that loess estimates are ordered by "x" for the plots, so that they are displayed properly
order.dist <- order(w.n$distance)
lines(x = w.n$distance[order.dist], y = predict(smooth.stand)[order.dist], lty = "solid", col = "red", lwd = 1)

# Standardized Pearson residual vs pi plot
plot(x = w.n$pi.hat, y = w.n$s.res, xlab = "Estimated probability of success", ylab = "Standardized Pearson residuals",
   main = "Standardized residuals vs. \n pi.hat", pch=chars)
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "blue")
smooth.stand <- loess(formula = s.res ~ pi.hat, data = w.n, weights = trials)
# Make sure that loess estimates are ordered by "X" for the plots, so that they are displayed properly
order.pi.hat <- order(w.n$pi.hat)
lines(x = w.n$pi.hat[order.pi.hat], y = predict(smooth.stand)[order.pi.hat], lty = "solid", col = "red", lwd = 1)

# Standardized Pearson residual vs linear predictor plot
plot(x = w.n$lin.pred, y = w.n$s.res, xlab = "Linear predictor", ylab = "Standardized Pearson residuals",
   main = "Standardized residuals vs. \n Linear predictor", pch=chars)
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "blue")
smooth.stand <- loess(formula = s.res ~ lin.pred, data = w.n, weights = trials)
# Make sure that loess estimates are ordered by "X" for the plots, so that they are displayed properly
order.lin.pred <- order(w.n$lin.pred)
lines(x = w.n$lin.pred[order.lin.pred], y = predict(smooth.stand)[order.lin.pred], lty = "solid", col = "red", lwd = 1)
# dev.off()  # Create plot for book


# Black-and-white version of plot
# pdf(file = "c:\\figures\\Figure5.3BW.pdf", width = 11, height = 6, colormodel = "cmyk", pointsize = 20)   # Create plot for book
par(mfrow = c(1,3))

# Standardized Pearson residual vs X plot
plot(x = w.n$distance, y = w.n$s.res, xlab = "Distance", ylab = "Standardized Pearson residuals",
   main = "Standardized residuals vs. \n X", pch=chars)
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "black")
lines(x = w.n$distance[order.dist], y = predict(smooth.stand)[order.dist], lty = "solid", col = "black", lwd = 1)

# Standardized Pearson residual vs pi plot
plot(x = w.n$pi.hat, y = w.n$s.res, xlab = "Estimated probability of success", ylab = "Standardized Pearson residuals",
   main = "Standardized residuals vs. \n pi.hat", pch=chars)
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "black")
lines(x = w.n$pi.hat[order.pi.hat], y = predict(smooth.stand)[order.pi.hat], lty = "solid", col = "black", lwd = 1)

# Standardized Pearson residual vs linear predictor plot
plot(x = w.n$lin.pred, y = w.n$s.res, xlab = "Linear predictor", ylab = "Standardized Pearson residuals",
   main = "Standardized residuals vs. \n Linear predictor", pch=chars)
abline(h = c(3, 2, 0, -2, -3), lty = "dotted", col = "black")
lines(x = w.n$lin.pred[order.lin.pred], y = predict(smooth.stand)[order.lin.pred], lty = "solid", col = "black", lwd = 1)
# dev.off()  # Create plot for book


##############################################################
# Goodness-of-Fit Measures and Tests
# 
# First the Deviance/DF
rdev <- mod.fit.bin$deviance 
dfr <- mod.fit.bin$df.residual 
ddf <- rdev/dfr 
thresh2 <- 1 + 2*sqrt(2/dfr) 
thresh3 <- 1 + 3*sqrt(2/dfr) 
c(rdev, dfr, ddf, thresh2, thresh3)
#
# Three goodness-of-fit tests are shown below: Hosmer and Lemeshow, Osius-Rojek, and Stukel. 
#  Each is contained in a separate function that assumes that a glm-class object 
#  has been created for the logistic regression. These functions are wrapped up into one 
#  function, AllGOFTests.R, that must be sourced.  All expect the model fit from glm to 
#  be in EVP form (there is no internal aggregation in the functions). Also, the 
#  OsiusRojek and Stukel expect that all interactions are actually listed as named 
#  individual variables (i.e. that they are included in the model as cross-product 
#  variables, like "X1:X2", rather than implicit interactions of other variables in the model.

# Before running, source the script:

source("AllGOFTests.R")  
HL <- HLTest(obj = mod.fit.bin, g = 10)
# Print out observed and expected counts in bins
cbind(HL$observed, round(HL$expect, digits = 1))
HL
# Pearson residuals for each group
round(HL$pear, digits = 1) 

o.r.test(obj = mod.fit.bin)

stukel.test(obj = mod.fit.bin)
