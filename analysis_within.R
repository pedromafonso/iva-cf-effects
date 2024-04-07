rm(list = ls())
# data =========================================================================
data <- readRDS("Data/data.rds")
data1 <- data$tw
data1$age_cat <- cut(data1$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                     labels = c("<12", "[12,18)", "≥18"))
remove(data)
# bias initiation date =========================================================
col <- rgb(0, 0, 0, alpha = 0.1)
title <- list(expression(alpha[0]~", intercept"),
              expression(alpha[1]~", time"),
              expression(gamma[0]~", change in intercept"),
              expression(gamma[1]~", change in slope"))
library("nlme")
m1a <- lme(GLI_FEV1_pct_predicted ~ time * ch_pt,
           random = ~ 1 | eDWID,
           correlation = corExp(form =~ time | eDWID, nugget = TRUE),
           data = data1)
int <- seq(0/12, 3/12, by = 0.5/12) # interval
#install.packages("foreach")
library("foreach")
#install.packages("doParallel")
library("doParallel")
n_cores <- parallel::detectCores() - 1
registerDoParallel(cores = n_cores)
est <- foreach(i = int, .packages= c("nlme")) %dopar% {
  bol <- abs(data1$time) >= i/2
  m <- try(update(m1a, data = data1[bol, ]), TRUE)
  if(!inherits(m, "try-error")) intervals(m, which = "fixed")[[1]] 
  else matrix(NA, p, 3) 
}
stopImplicitCluster()
est <- array(unlist(est), dim = c(nrow(est[[1]]), ncol(est[[1]]), length(est)),
             dimnames = list(rownames(est[[1]]), colnames(est[[1]]), int * 12))
png("manuscript/figures/initiation_bias.png", width = 10 * 300, height = 6 * 300, res = 300)
{
  par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1), bg = "white")
  for(i in 1:4) {
    ylim <- range(est[i, , ], na.rm = TRUE)
    xlim <- range(int * 12)
    plot(NA, type = 'n', xlim = xlim, ylim = ylim, 
         bty = "n", xaxt = "n", yaxt = "n", 
         ylab = "Estimate", xlab = "Exclusion interval length (months)")
    box(lwd = 0.5)
    axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
    axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
    # polygon & line
    x <- c(int * 12, rev(int * 12))
    y <- c(est[i, 1, ], rev(est[i, 3, ]))
    polygon(x = x[!is.na(y)], y = y[!is.na(y)], col = col, border = NA, lwd = 1)
    lines(int * 12, est[i, 2, ], lwd = 1, type = "b", pch = 20)
    abline(h = 0, lwd = 0.5)
    title(title[[i]], line = 1)
  }
}
dev.off()
n0 <- length(unique(data1$eDWID)); n0 # 560
N0 <- nrow(data1); N0 # 20424
bol <- abs(data1$time) >= (0.5/12)/2
data1 <- data1[bol, ]
n1 <- length(unique(data1$eDWID)); n1 # 560
N1 <- nrow(data1); N1 # 19860
N0 - N1 # 564
(N0 - N1)/N0 * 100 #  2.761457
saveRDS(data1, "Data/data1.rds")
# tables =======================================================================
quantile2 <- function(x, digits = 2) {
  qts <- quantile(x, probs = c(0.25, 0.5, 0.75))[c(2, 1, 3)]
  qts <- format(round(qts, digits), nsmall = digits)
  qts <- gsub(" ", "", qts, fixed = TRUE) # remove empty spaces from the qts output
  cat(paste0(qts[1], " (", qts[2], "–", qts[3], ")", collapse = ""))
}
table2 <- function(x, digits = 2) {
  tb <- table(x, exclude = NULL)
  pc <- format(round(tb/length(x) * 100, digits), nsmall = digits)
  pc <- gsub(" ", "", pc, fixed = TRUE) # remove empty spaces from the pc output
  res <- paste0(tb, " (", pc, "%)")
  names(res) <- names(tb)
  cat(res)
}
## table main - treatment between ==============================================
### number individuals
length(unique(data1$eDWID)) # [1] 560
### number measurements
nrow(data1) # [1] 20424
### number measurements per individual
quantile2(table(data1$eDWID)) # 32.50 (21.00–46.25)
### vising times
png("manuscript/figures/vtimes_trt_w.png", width = 5 * 300, height = 4 * 300, res = 300) 
{
  par(bg = "white")
  hist(data1$time, main = "", freq = FALSE, 
       xlab = "Visiting time (years)", bty = "n", xaxt = "n", yaxt = "n", 
       col = rgb(0, 0, 0, 0.5), border = NA)
  box(lwd = 0.5)
  axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
  axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
}
dev.off()
### total follow-up duration
head_rows <- tapply(seq_along(data1$eDWID), data1$eDWID, head, 1)
tail_rows <- tapply(seq_along(data1$eDWID), data1$eDWID, tail, 1)
sum(data1$encounterage[tail_rows] - data1$encounterage[head_rows]) # [1] 4449.937
remove(head_rows, tail_rows)
### follow-up duration per individual
head_rows <- tapply(seq_along(data1$eDWID), data1$eDWID, head, 1)
tail_rows <- tapply(seq_along(data1$eDWID), data1$eDWID, tail, 1)
quantile2(data1$encounterage[tail_rows] - data1$encounterage[head_rows]) # 9.05 (5.89–10.24)
remove(head_rows, tail_rows)
# within-subject analysis ======================================================
col <- rgb(0, 0, 0, alpha = 0.1)
library(nlme)
m1a <- lme(GLI_FEV1_pct_predicted ~ (time + I(time * (time > 0)) + ch_pt) * age_cat,
           random = ~ 1 | eDWID,
           correlation = corExp(form =~ time | eDWID, nugget = TRUE),
           data = data1, method = "ML")
library(splines)
# m1a <- lme(GLI_FEV1_pct_predicted ~ (ns(time, df = 1, B = c(-6, 0))
#            + I(ns(time, df = 1, B = c(0, 6)) * (time > 0)) + ch_pt) * age_cat,
#            random = ~ 1 | eDWID,
#            correlation = corExp(form =~ time | eDWID, nugget = TRUE),
#            data = data1, method = "ML")
m1b <- lme(GLI_FEV1_pct_predicted ~ (ns(time, df = 2, B = c(-6, 0)) 
           + I(ns(time, df = 2, B = c(0, 6)) * (time > 0)) + ch_pt) * age_cat,
           random = ~ 1 | eDWID,
           correlation = corExp(form =~ time | eDWID, nugget = TRUE),
           data = data1, method = "ML")
m1c <- lme(GLI_FEV1_pct_predicted ~ (ns(time, df = 3, B = c(-6, 0)) 
           + I(ns(time, df = 3, B = c(0, 6)) * (time > 0)) + ch_pt) * age_cat,
           random = ~ 1 | eDWID,
           correlation = corExp(form =~ time | eDWID, nugget = TRUE),
           data = data1, method = "ML")
## table - model selection =====================================================
AIC(m1a, m1b, m1c)
BIC(m1a, m1b, m1c)
anova(m1a, m1b) # 0.0337
anova(m1a, m1c) # 0.1022
## plot - model selection ======================================================
get_est <- function(object, res = 100) {
  toy <- expand.grid(time = seq(-6, 6, length.out = res),
                     age_cat = factor(c("<12", "[12,18)", "≥18"), 
                                      levels = c("<12", "[12,18)", "≥18")))
  toy$ch_pt <- as.numeric(toy$time > 0)
  frml <- formula(object)[c(1, 3)]
  X <- model.matrix(frml, data = toy)
  betas <- if(inherits(object, "lme")) {
    nlme::fixef(object)
  } else if(inherits(object, "gee")) {
    object$coefficients
  }
  D <- vcov(object)
  se <- apply(X, 1, function(x) sqrt(t(x) %*% D %*% x))
  toy$mean <- X %*% betas # mean
  toy$upper <- toy$mean + qnorm(0.975) * se # upper
  toy$lower <- toy$mean - qnorm(0.975) * se # lower
  toy
}
est1a <- get_est(m1a)
est1b <- get_est(m1b)
est1c <- get_est(m1c)
png("manuscript/figures/model_selection1.png", width = 6 * 300, height = 8 * 300, res = 300)
{
  ylim <- range(est1a$lower, est1a$upper, est1b$lower, est1b$upper, est1c$lower, 
                est1c$upper)
  par(mfrow = c(3, 1), bg = "white", mar = c(5.1, 4.1, 4.1 - 2, 2.1))
  for(i in c("<12", "[12,18)", "≥18")) {
    bol <- est1a$age_cat == i
    plot(NA, xlim = c(-6, 6), ylim = ylim, type = 'n', bty = "n", 
         xaxt = "n", yaxt = "n", 
         ylab = expression(ppFEV[1]), 
         xlab = "Time since ivacaftor initiation (years)")
    mtext(i, side = 3, line = 0.5, cex = 0.8)
    box(lwd = 0.5)
    axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
    axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
    bol1 <- bol & est1a$time <= 0
    bol2 <- bol & est1a$time >= 0
    polygon(x = c(est1a$time[bol1], rev(est1a$time[bol1])), 
            y = c(est1a$lower[bol1], rev(est1a$upper[bol1])), col = col, border = NA)
    polygon(x = c(est1a$time[bol2], rev(est1a$time[bol2])), 
            y = c(est1a$lower[bol2], rev(est1a$upper[bol2])), col = col, border = NA)
    polygon(x = c(est1b$time[bol1], rev(est1b$time[bol1])), 
            y = c(est1b$lower[bol1], rev(est1b$upper[bol1])), col = col, border = NA)
    polygon(x = c(est1b$time[bol2], rev(est1b$time[bol2])), 
            y = c(est1b$lower[bol2], rev(est1b$upper[bol2])), col = col, border = NA)
    polygon(x = c(est1c$time[bol1], rev(est1c$time[bol1])), 
            y = c(est1c$lower[bol1], rev(est1c$upper[bol1])), col = col, border = NA)
    polygon(x = c(est1c$time[bol2], rev(est1c$time[bol2])), 
            y = c(est1c$lower[bol2], rev(est1c$upper[bol2])), col = col, border = NA)
    lines(est1a$time[bol1], est1a$mean[bol1])
    lines(est1a$time[bol2], est1a$mean[bol2])
    lines(est1b$time[bol1], est1b$mean[bol1])
    lines(est1b$time[bol2], est1b$mean[bol2])
    lines(est1c$time[bol1], est1c$mean[bol1])
    lines(est1c$time[bol2], est1c$mean[bol2])
    bol <- data1$age_cat == i
    segments(x0 = data1$time[bol], x1 = data1$time[bol], y0 = par("usr")[3], 
             y1 = par("usr")[3] + 0.05 * diff(par("usr")[3:4]), col = col)
  }
}
dev.off()
m1a <- update(m1a, method = "REML")
m1d <- lme(GLI_FEV1_pct_predicted ~ (time + I(time * (time > 0)) + ch_pt) * age_cat,
           random = ~ 1 | eDWID/ch_pt,
           correlation = corExp(form =~ time | eDWID/ch_pt, nugget = TRUE),
           data = data1)
p <- 0.5 * (1 - pchisq(anova(m1a, m1d)$L.Ratio[2], 1)); p
m1 <- m1a
m1 <- update(m1, fixed = GLI_FEV1_pct_predicted ~ -1 + time:ch_pt:age_cat + age_cat:ch_pt)
round(intervals(m1, which = "fixed")[[1]][, c(2, 1, 3)], 2) # estimates before vs after
m1 <- update(m1, fixed = GLI_FEV1_pct_predicted ~ -1 + (time * ch_pt) : age_cat + age_cat)
round(intervals(m1, which = "fixed")[[1]][, c(2, 1, 3)], 2) # ivacaftor effect
## effects plot ================================================================
png("manuscript/figures/effects1.png", width = 8 * 300, height = 5 * 300, res = 300)
{
  res <- 3
  toy1 <- expand.grid(time = seq(-6, 0, length.out = res),
                      age_cat = factor(c("<12", "[12,18)", "≥18"), 
                                       levels = c("<12", "[12,18)", "≥18")),
                      ch_pt = factor(c("1"), levels = c("1", "2")))
  toy2 <- expand.grid(time = seq(0, 6, length.out = res),
                      age_cat = factor(c("<12", "[12,18)", "≥18"), 
                                       levels = c("<12", "[12,18)", "≥18")),
                      ch_pt = factor(c("1", "2"), levels = c("1", "2")))
  frml <- formula(m1)[c(1, 3)]
  X1 <- model.matrix(frml, data = toy1)
  X2 <- model.matrix(frml, data = toy2)
  betas <- fixef(m1)
  D <- vcov(m1)
  se1 <- apply(X1, 1, function(x) sqrt(t(x) %*% D %*% x))
  se2 <- apply(X2, 1, function(x) sqrt(t(x) %*% D %*% x))
  toy1$mean <- X1 %*% betas # mean
  toy1$upper <- toy1$mean + qnorm(0.975) * se1 # upper
  toy1$lower <- toy1$mean - qnorm(0.975) * se1 # lower
  toy2$mean <- X2 %*% betas # mean
  toy2$upper <- toy2$mean + qnorm(0.975) * se2 # upper
  toy2$lower <- toy2$mean - qnorm(0.975) * se2 # lower
  ylim <- range(toy1$lower, toy1$upper, toy2$lower, toy2$upper)
  par(mfrow = c(1, 3), bg = "white", mar = c(5.1, 4.1, 4.1 - 2, 2.1))
  for(i in c("<12", "[12,18)", "≥18")) {
    plot(NA, xlim = c(-6, 6), ylim = ylim, type = 'n', bty = "n", 
         xaxt = "n", yaxt = "n", 
         ylab = expression(ppFEV[1]), 
         xlab = "Time since ivacaftor initiation (years)")
    mtext(i, side = 3, line = 0.5, cex = 0.8)
    box(lwd = 0.5)
    abline(v = 0, lwd = 0.5)
    axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
    axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
    bol1 <- toy1$age_cat == i
    bol21 <- toy2$age_cat == i & toy2$ch_pt == "1"
    bol22 <- toy2$age_cat == i & toy2$ch_pt == "2"
    polygon(x = c(toy1$time[bol1], rev(toy1$time[bol1])), 
            y = c(toy1$lower[bol1], rev(toy1$upper[bol1])), col = col, border = NA)
    polygon(x = c(toy2$time[bol21], rev(toy2$time[bol21])), 
            y = c(toy2$lower[bol21], rev(toy2$upper[bol21])), col = col, border = NA)
    polygon(x = c(toy2$time[bol22], rev(toy2$time[bol22])), 
            y = c(toy2$lower[bol22], rev(toy2$upper[bol22])), col = col, border = NA)
    lines(toy1$time[bol1], toy1$mean[bol1])
    lines(toy2$time[bol21], toy2$mean[bol21], lty = 4)
    lines(toy2$time[bol22], toy2$mean[bol22])
    legend("bottomleft", legend = c("Natural progression"),
           lwd = 1, lty = c(4), bty = "n")
  }
}
dev.off()
# png("manuscript/figures/effects1.png", width = 8 * 300, height = 5 * 300, res = 300)
# {
#   res <- 3
#   toy <- expand.grid(time = seq(0, 6, length.out = res),
#                      age_cat = factor(c("<12", "[12,18)", "≥18"), 
#                                       levels = c("<12", "[12,18)", "≥18")),
#                      ch_pt = factor(c("1", "2"), levels = c("1", "2")))
#   frml <- formula(m1)[c(1, 3)]
#   X <- model.matrix(frml, data = toy)
#   betas <- fixef(m1)
#   D <- vcov(m1)
#   se <- apply(X, 1, function(x) sqrt(t(x) %*% D %*% x))
#   toy$mean <- X %*% betas # mean
#   toy$upper <- toy$mean + qnorm(0.975) * se # upper
#   toy$lower <- toy$mean - qnorm(0.975) * se # lower
#   ylim <- range(toy$lower, toy$upper)
#   par(mfrow = c(1, 3), bg = "white", mar = c(5.1, 4.1, 4.1 - 2, 2.1))
#   for(i in c("<12", "[12,18)", "≥18")) {
#     plot(NA, xlim = c(0, 6), ylim = ylim, type = 'n', bty = "n", 
#          xaxt = "n", yaxt = "n", 
#          ylab = expression(ppFEV[1]), 
#          xlab = "Time since ivacaftor initiation (years)")
#     mtext(i, side = 3, line = 0.5, cex = 0.8)
#     box(lwd = 0.5)
#     axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
#     axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
#     bol1 <- toy$age_cat == i & toy$ch_pt == "1"
#     bol2 <- toy$age_cat == i & toy$ch_pt == "2"
#     polygon(x = c(toy$time[bol1], rev(toy$time[bol1])), 
#             y = c(toy$lower[bol1], rev(toy$upper[bol1])), col = col, border = NA)
#     polygon(x = c(toy$time[bol2], rev(toy$time[bol2])), 
#             y = c(toy$lower[bol2], rev(toy$upper[bol2])), col = col, border = NA)
#     lines(toy$time[bol1], toy$mean[bol1], lty = 4)
#     lines(toy$time[bol2], toy$mean[bol2])
#     legend("bottomleft", legend = c("Natural progression", "On ivacaftor"),
#            lwd = 1, lty = c(4, 1), bty = "n")
#   }
# }
# dev.off()
## residuals ===================================================================
col <- rgb(0, 0, 0, alpha = 0.1)
res <- residuals(m1, type = "pearson")
fit <- fitted(m1)
png("manuscript/figures/residuals1.png", width = 10 * 300, height = 6 * 300, res = 300)
{
  par(mfrow = c(1, 2), bg = "white")
  plot(NA, xlim = range(fit), ylim = range(res), bty = "n", xaxt = "n",  yaxt = "n",
       ylab = "Standardized residual", xlab = "Fitted value")
  axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
  axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
  box(lwd = 0.5)
  abline(h = 0, lwd = 0.5)
  points(fit, res, pch = 21, col = NA, bg = col)
  qqnorm(res, bty = "n", xaxt = "n", yaxt = "n",
         pch = 21, col = NA, bg = col, main = "",
         ylab = "Quantiles of standard normal", xlab = "Standardized residual")
  qqline(res, lwd = 0.5)
  box(lwd = 0.5)
  axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
  axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
}
dev.off()
