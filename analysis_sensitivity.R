rm(list = ls())
# data =========================================================================
data1 <- readRDS("Data/data1.rds")
data2 <- readRDS("Data/data2.rds")
# model ========================================================================
## within analysis =============================================================
library(nlme)
m1 <- lme(GLI_FEV1_pct_predicted ~ -1 + time:ch_pt:age_cat + age_cat:ch_pt,
          random = ~ 1 | eDWID,
          correlation = corExp(form =~ time | eDWID, nugget = TRUE),
          data = data1)
m1_3y <- update(m1, data = data1[abs(data1$time) <= 3, ])
round(intervals(m1_3y, which = "fixed")[[1]][, c(2, 1, 3)], 2) # estimates before vs after
m1_3y <- update(m1_3y, fixed = GLI_FEV1_pct_predicted ~ -1 + (time * ch_pt) : age_cat + age_cat)
round(intervals(m1_3y, which = "fixed")[[1]][, c(2, 1, 3)], 2) # ivacaftor effect
## between analysis ============================================================
m2 <- lme(GLI_FEV1_pct_predicted ~ -1 + time:group:age_cat + age_cat:group,
          random = ~ 1 | subclass/eDWID,
          correlation = corExp(form =~ time | subclass/eDWID, nugget = TRUE),
          data = data2)
m2_3y <- update(m2, data = data2[abs(data2$time) <= 3, ], 
                control = list(msMaxIter = 1000, msMaxEval = 1000))
round(intervals(m2_3y, which = "fixed")[[1]][, c(2, 1, 3)], 2) # estimates control vs treatment
m2_3y <- update(m2_3y, fixed = GLI_FEV1_pct_predicted ~ -1 + (time * group) : age_cat + age_cat)
round(intervals(m2_3y, which = "fixed")[[1]][, c(2, 1, 3)], 2) # ivacaftor effect
# effects plot =================================================================
col <- rgb(0, 0, 0, alpha = 0.1)
## within analysis =============================================================
png("manuscript/figures/effects1_3y.png", width = 8 * 300, height = 5 * 300, res = 300)
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
  betas_3y <- fixef(m1_3y)
  D <- vcov(m1)
  D_3y <- vcov(m1_3y)
  se1 <- apply(X1, 1, function(x) sqrt(t(x) %*% D %*% x))
  se2 <- apply(X2, 1, function(x) sqrt(t(x) %*% D %*% x))
  se1_3y <- apply(X1, 1, function(x) sqrt(t(x) %*% D_3y %*% x))
  se2_3y <- apply(X2, 1, function(x) sqrt(t(x) %*% D_3y %*% x))
  toy1$mean <- X1 %*% betas # mean
  toy1$mean_3y <- X1 %*% betas_3y # mean
  toy1$upper <- toy1$mean + qnorm(0.975) * se1 # upper
  toy1$upper_3y <- toy1$mean_3y + qnorm(0.975) * se1_3y # upper
  toy1$lower <- toy1$mean - qnorm(0.975) * se1 # lower
  toy1$lower_3y <- toy1$mean_3y - qnorm(0.975) * se1_3y # lower
  toy2$mean <- X2 %*% betas # mean
  toy2$mean_3y <- X2 %*% betas_3y # mean
  toy2$upper <- toy2$mean + qnorm(0.975) * se2 # upper
  toy2$upper_3y <- toy2$mean_3y + qnorm(0.975) * se2_3y # upper
  toy2$lower <- toy2$mean - qnorm(0.975) * se2 # lower
  toy2$lower_3y <- toy2$mean_3y - qnorm(0.975) * se2_3y # lower
  ylim <- range(toy1$lower, toy1$upper, toy2$lower, toy2$upper,
                toy1$lower_3y, toy1$upper_3y, toy2$lower_3y, toy2$upper_3y)
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
    bol1_3y <- toy1$age_cat == i & toy1$time >= -3
    bol21 <- toy2$age_cat == i & toy2$ch_pt == "1"
    bol21_3y <- toy2$age_cat == i & toy2$ch_pt == "1" & toy2$time <= 3
    bol22 <- toy2$age_cat == i & toy2$ch_pt == "2"
    bol22_3y <- toy2$age_cat == i & toy2$ch_pt == "2" & toy2$time <= 3
    polygon(x = c(toy1$time[bol1], rev(toy1$time[bol1])), 
            y = c(toy1$lower[bol1], rev(toy1$upper[bol1])), col = col, border = NA)
    polygon(x = c(toy1$time[bol1_3y], rev(toy1$time[bol1_3y])), 
            y = c(toy1$lower_3y[bol1_3y], rev(toy1$upper_3y[bol1_3y])), col = col, border = NA)
    polygon(x = c(toy2$time[bol21], rev(toy2$time[bol21])), 
            y = c(toy2$lower[bol21], rev(toy2$upper[bol21])), col = col, border = NA)
    polygon(x = c(toy2$time[bol21_3y], rev(toy2$time[bol21_3y])), 
            y = c(toy2$lower_3y[bol21_3y], rev(toy2$upper_3y[bol21_3y])), col = col, border = NA)
    
    polygon(x = c(toy2$time[bol22], rev(toy2$time[bol22])), 
            y = c(toy2$lower[bol22], rev(toy2$upper[bol22])), col = col, border = NA)
    polygon(x = c(toy2$time[bol22_3y], rev(toy2$time[bol22_3y])), 
            y = c(toy2$lower_3y[bol22_3y], rev(toy2$upper_3y[bol22_3y])), col = col, border = NA)
    lines(toy1$time[bol1], toy1$mean[bol1])
    lines(toy1$time[bol1_3y], toy1$mean_3y[bol1_3y])
    lines(toy2$time[bol21], toy2$mean[bol21], lty = 4)
    lines(toy2$time[bol21_3y], toy2$mean_3y[bol21_3y], lty = 4)
    lines(toy2$time[bol22], toy2$mean[bol22])
    lines(toy2$time[bol22_3y], toy2$mean_3y[bol22_3y])
    legend("bottomleft", legend = c("Natural progression"),
           lwd = 1, lty = c(4), bty = "n")
  }
}
dev.off()
## between analysis ============================================================
png("manuscript/figures/effects2_3y.png", width = 8 * 300, height = 5 * 300, res = 300)
{
  res <- 3
  toy <- expand.grid(time = seq(0, 6, length.out = res),
                     age_cat = factor(c("<12", "[12,18)", "≥18"), 
                                      levels = c("<12", "[12,18)", "≥18")),
                     group = factor(c("ctr", "trt"), levels = c("ctr", "trt")))
  frml <- formula(m2)[c(1, 3)]
  X <- model.matrix(frml, data = toy)
  betas <- fixef(m2)
  betas_3y <- fixef(m2_3y)
  D <- vcov(m2)
  D_3y <- vcov(m2_3y)
  se <- apply(X, 1, function(x) sqrt(t(x) %*% D %*% x))
  se_3y <- apply(X, 1, function(x) sqrt(t(x) %*% D_3y %*% x))
  toy$mean <- X %*% betas # mean
  toy$mean_3y <- X %*% betas_3y # mean
  toy$upper <- toy$mean + qnorm(0.975) * se # upper
  toy$upper_3y <- toy$mean_3y + qnorm(0.975) * se_3y # upper
  toy$lower <- toy$mean - qnorm(0.975) * se # lower
  toy$lower_3y <- toy$mean_3y - qnorm(0.975) * se_3y # lower
  ylim <- range(toy$lower, toy$upper, toy$lower_3y, toy$upper_3y)
  par(mfrow = c(1, 3), bg = "white", mar = c(5.1, 4.1, 4.1 - 2, 2.1))
  for(i in c("<12", "[12,18)", "≥18")) {
    plot(NA, xlim = c(0, 6), ylim = ylim, type = 'n', bty = "n", 
         xaxt = "n", yaxt = "n", 
         ylab = expression(ppFEV[1]), 
         xlab = "Time since ivacaftor initiation (years)")
    mtext(i, side = 3, line = 0.5, cex = 0.8)
    box(lwd = 0.5)
    axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
    axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
    bol1 <- toy$age_cat == i & toy$group == "ctr"
    bol1_3y <- toy$age_cat == i & toy$group == "ctr" & toy$time <= 3
    bol2 <- toy$age_cat == i & toy$group == "trt"
    bol2_3y <- toy$age_cat == i & toy$group == "trt" & toy$time <= 3
    polygon(x = c(toy$time[bol1], rev(toy$time[bol1])), 
            y = c(toy$lower[bol1], rev(toy$upper[bol1])), col = col, border = NA)
    polygon(x = c(toy$time[bol1_3y], rev(toy$time[bol1_3y])), 
            y = c(toy$lower_3y[bol1_3y], rev(toy$upper_3y[bol1_3y])), col = col, border = NA)
    polygon(x = c(toy$time[bol2], rev(toy$time[bol2])), 
            y = c(toy$lower[bol2], rev(toy$upper[bol2])), col = col, border = NA)
    polygon(x = c(toy$time[bol2_3y], rev(toy$time[bol2_3y])), 
            y = c(toy$lower_3y[bol2_3y], rev(toy$upper_3y[bol2_3y])), col = col, border = NA)
    lines(toy$time[bol1], toy$mean[bol1], lty = 2)
    lines(toy$time[bol1_3y], toy$mean_3y[bol1_3y], lty = 2)
    lines(toy$time[bol2], toy$mean[bol2])
    lines(toy$time[bol2_3y], toy$mean_3y[bol2_3y])
    legend("bottomleft", legend = c("Control", "On ivacaftor"),
           lwd = 1, lty = 2:1, bty = "n")
  }
}
dev.off()