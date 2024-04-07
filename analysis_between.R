rm(list = ls())
# data =========================================================================
data <- readRDS("Data/data.rds")
trt0 <- data$tb0
ctr0 <- data$c0
trt  <- data$tb
ctr  <- data$c
remove(data)
trt0$group <- factor("trt", levels = c("ctr", "trt"))
ctr0$group <- factor("ctr", levels = c("ctr", "trt"))
match_vars <- c("ann_aztreonam", "cfrd_status2", "ann_colistin", 
                "ann_corticosteroids1", "ann_corticosteroids2", "Diagnosis_year",
                "ann_dornasealfa", "encounterage", "ann_fungalyeast1", "Gender",
                "ann_GLI_FEV1_pct_predicted", "ann_GLI_FVC_pct_predicted",
                "ann_GLI_FEV1FVC_pct_predicted", "ann_leukotriene_others3",
                "ann_mycobacterials_nontub", "ann_pseudomonasaeruginosa",
                "ReviewYear", "ann_staphyl_type1", "ann_tobi")
match_data <- rbind(ctr0[, c("group", "eDWID", match_vars)], 
                    trt0[, c("group", "eDWID", match_vars)])
data2 <- rbind(ctr[, !names(ctr) %in% "stable2"], 
               trt[, !names(trt) %in% "iva_st"])

bol_ctr <- data2$eDWID %in% unique(match_data$eDWID[match_data$group == "ctr"]) 
bol_trt <- data2$eDWID %in% unique(match_data$eDWID[match_data$group == "trt"])
n0_ctr <- length(unique(data2$eDWID[bol_ctr])); n0_ctr # 10294
n0_trt <- length(unique(data2$eDWID[bol_trt])); n0_trt # 560
N0_ctr <- nrow(data2[bol_ctr, ]); N0_ctr # 231678
N0_trt <- nrow(data2[bol_trt, ]); N0_trt # 9233
bol <- abs(data2$time) >= (0.5/12)/2
data2 <- data2[bol, ]
bol_ctr <- data2$eDWID %in% unique(match_data$eDWID[match_data$group == "ctr"]) 
bol_trt <- data2$eDWID %in% unique(match_data$eDWID[match_data$group == "trt"])
n1_ctr <- length(unique(data2$eDWID[bol_ctr])); n1_ctr # 10294
n1_trt <- length(unique(data2$eDWID[bol_trt])); n1_trt # 560
N1_ctr <- nrow(data2[bol_ctr, ]); N1_ctr # 231588
N1_trt <- nrow(data2[bol_trt, ]); N1_trt # 9232
N0_ctr - N1_ctr + n1_ctr # 10384 # n1_ctr to account for all zeros we excluded before, one per individual
N0_trt - N1_trt + n1_trt # 561 # n1_trt to account for all zeros we excluded before, one per individual
(N0_ctr - N1_ctr - n1_ctr)/N0_ctr * 100 #  4.404389
(N0_trt - N1_trt - n1_trt)/N0_trt * 100 #  6.05437
# matching =====================================================================
match_data$group2 <- as.numeric(match_data$group == "trt")
match_form <- as.formula(paste("group ~ ", paste0(match_vars, collapse = "+")))
#install.packages("MatchIt")
library(MatchIt)
match_obj <- matchit(formula = match_form, data = match_data, 
                     method = "nearest", distance = "glm", link = 'logit', 
                     replace = FALSE, ratio = 5)
match_data2 <- match.data(match_obj)
data2 <- data2[data2$eDWID %in% unique(match_data2$eDWID), ]
ctr20 <- ctr0[ctr0$eDWID %in% unique(match_data2$eDWID), ]
ctr2  <- ctr[ctr$eDWID %in% unique(match_data2$eDWID), ]
data2 <- merge(data2, match_data2[, c("eDWID", "subclass", "group")], 
               by = "eDWID", all.x = TRUE)
data2$age_cat <- cut(data2$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                     labels = c("<12", "[12,18)", "≥18"))
remove(match_form, match_vars)
## covariate balance plot ======================================================
new_names <- c(ann_aztreonam = "Aztreonam", 
               cfrd_status2_0 = "CF-related diabetes (no)", 
               cfrd_status2_2 = "CF-related diabetes (impaired)", 
               cfrd_status2_345 = "CF-related diabetes (type I–II)", 
               ann_colistin = "Colistin", 
               ann_corticosteroids1 = "Corticosteroids (oral)", 
               ann_corticosteroids2 = "Corticosteroids (inhaled)", 
               Diagnosis_year = "Diagnosis year",
               ann_dornasealfa = "Dornase alfa", 
               encounterage = "Age", 
               ann_fungalyeast1 = "Aspergillus fumigatus", 
               Gender_M = "Sex (male)",
               ann_GLI_FEV1_pct_predicted = "ppFEV1", 
               ann_GLI_FVC_pct_predicted = "ppFVC",
               ann_GLI_FEV1FVC_pct_predicted = "ppFEV1/FVC", 
               ann_leukotriene_others3 = "Antifungals",
               ann_mycobacterials_nontub = "Nontuberculous mycobacterium", 
               ann_pseudomonasaeruginosa = "Pseudomonas aeruginosa",
               ReviewYear = "Encounter year", 
               ann_staphyl_type1 = "Methicillin-resistant Staphylococcus aureus", 
               ann_tobi = "Tobramycin")
#install.packages("cobalt")
library(cobalt)
#install.packages(ggplot2)
library(ggplot2)
png("manuscript/figures/matching.png", width = 500 * 6, height = 300 * 6, res = 300)
{
  love.plot(match_obj, binary = "std", title = "",
            shapes = c("circle filled", "circle"), # pch
            drop.distance = TRUE,
            sample.names = c("Unmatched", "PS matched"),
            var.names = new_names,
            thresholds = c(m = .1)) + 
    xlab("Standardized mean differences") +
    theme(legend.title = element_blank())
}
dev.off()
# tables =======================================================================
saveRDS(data2, "Data/data2.rds")
quantile2 <- function(x, digits = 2) {
  qts <- quantile(x, probs = c(0.25, 0.5, 0.75))[c(2, 1, 3)]
  qts <- format(round(qts, digits), nsmall = digits)
  qts <- gsub(" ", "", qts, fixed = TRUE) # remove empty spaces from the qts output
  cat(paste0(qts[1], " (", qts[2], "–", qts[3], ")", collapse = ""))
}
table2 <- function(x, digits = 2) {
  tb <- table(x, exclude = NULL)
  pc <- format(round(tb/length(x)*100, digits), nsmall = digits)
  pc <- gsub(" ", "", pc, fixed = TRUE) # remove empty spaces from the pc output
  res <- paste0(tb, " (", pc, "%)")
  names(res) <- names(tb)
  cat(res)
}
## table main - control ========================================================
### number individuals
length(unique(data2$eDWID[data2$group == "ctr"])) # [1] 2800
### number measurements
nrow(data2[data2$group == "ctr", ]) # [1] 48913
quantile2(table(data2$eDWID[data2$group == "ctr"])) # 15.00 (9.00–23.00)
### vising times
png("manuscript/figures/vtimes_ctr.png", width = 5 * 300, height = 4 * 300, res = 300) 
{
  par(bg = "white")
  hist(data2$time[data2$group == "ctr"], main = "", freq = FALSE, 
       xlab = "Visiting time (years)", bty = "n", xaxt = "n", yaxt = "n", 
       col = rgb(0, 0, 0, 0.5), border = NA)
  box(lwd = 0.5)
  axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
  axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
}
dev.off()
### total follow-up duration
head_rows <- tapply(seq_along(data2$eDWID[data2$group == "ctr"]), 
                    data2$eDWID[data2$group == "ctr"], head, 1)
tail_rows <- tapply(seq_along(data2$eDWID[data2$group == "ctr"]), 
                    data2$eDWID[data2$group == "ctr"], tail, 1)
sum(data2$encounterage[data2$group == "ctr"][tail_rows] - 
      data2$encounterage[data2$group == "ctr"][head_rows]) # [1] 8862.587
remove(head_rows, tail_rows)
### follow-up duration per individual
head_rows <- tapply(seq_along(data2$eDWID[data2$group == "ctr"]), 
                    data2$eDWID[data2$group == "ctr"], head, 1)
tail_rows <- tapply(seq_along(data2$eDWID[data2$group == "ctr"]), 
                    data2$eDWID[data2$group == "ctr"], tail, 1)
quantile2(data2$encounterage[data2$group == "ctr"][tail_rows] - 
            data2$encounterage[data2$group == "ctr"][head_rows]) # 3.18 (2.23–4.14)
remove(head_rows, tail_rows)
### matching encounter year
quantile2(ctr20$ReviewYear) # 2012 (2012–2013)
### ppFEV1 matching encounter 
quantile2(ctr20$ann_GLI_FEV1_pct_predicted) # 87.28 (68.91–100.75)
### age matching encounter
quantile2(ctr20$encounterage) # 15.10 (7.01–24.78)
### age_cat matching encounter
ctr20$age_cat <- cut(ctr20$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                     labels = c("<12", "[12,18)", "≥18"))
table2(ctr20$age_cat) # 1166 (41.64%) 466 (16.64%) 1168 (41.71%)
### lung transplant
table(data2$txp_flag21[data2$group == "ctr"]) # 52
round(table(data2$txp_flag21[data2$group == "ctr"])/2800 * 100, 2) # 1.86
### death 
table(data2$death_flag21[data2$group == "ctr"]) # 13
round(table(data2$death_flag21[data2$group == "ctr"])/2800 * 100, 2) # 0.46
## table main - treatment between ==============================================
### number individuals
length(unique(data2$eDWID[data2$group == "trt"])) # [1] 560
### number measurements
nrow(data2[data2$group == "trt", ]) # [1] 9232
quantile2(table(data2$eDWID[data2$group == "trt"])) # 15.00 (10.00–20.00)
### vising times
png("manuscript/figures/vtimes_trt_b.png", width = 5 * 300, height = 4 * 300, res = 300) 
{
  par(bg = "white")
  hist(data2$time[data2$group == "trt"], main = "", freq = FALSE, 
       xlab = "Visiting time (years)", bty = "n", xaxt = "n", yaxt = "n", 
       col = rgb(0, 0, 0, 0.5), border = NA)
  box(lwd = 0.5)
  axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
  axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
}
dev.off()
### total follow-up duration
head_rows <- tapply(seq_along(data2$eDWID[data2$group == "trt"]), 
                    data2$eDWID[data2$group == "trt"], head, 1)
tail_rows <- tapply(seq_along(data2$eDWID[data2$group == "trt"]), 
                    data2$eDWID[data2$group == "trt"], tail, 1)
sum(data2$encounterage[data2$group == "trt"][tail_rows] - 
      data2$encounterage[data2$group == "trt"][head_rows]) # [1] 2048.538
remove(head_rows, tail_rows)
### follow-up duration per individual
head_rows <- tapply(seq_along(data2$eDWID[data2$group == "trt"]), 
                    data2$eDWID[data2$group == "trt"], head, 1)
tail_rows <- tapply(seq_along(data2$eDWID[data2$group == "trt"]), 
                    data2$eDWID[data2$group == "trt"], tail, 1)
quantile2(data2$encounterage[data2$group == "trt"][tail_rows] - 
            data2$encounterage[data2$group == "trt"][head_rows]) # 4.04 (3.20–4.47)
remove(head_rows, tail_rows)
### matching encounter year
quantile2(trt0$ReviewYear) # 2012.00 (2012.00–2012.00)
### ppFEV1 matching encounter 
quantile2(trt0$ann_GLI_FEV1_pct_predicted) # 88.44 (70.60–99.99)
### age matching encounter
quantile2(trt0$encounterage) # 15.62 (9.06–23.98)
### age_cat matching encounter
trt0$age_cat <- cut(trt0$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                    labels = c("<12", "[12,18)", "≥18"))
table2(trt0$age_cat) # 199 (35.54%) 133 (23.75%) 228 (40.71%)
### lung transplant
table(data2$txp_flag21[data2$group == "trt"]) # 2
round(table(data2$txp_flag21[data2$group == "trt"])/560 * 100, 2) # 0.36
### death 
table(data2$death_flag21[data2$group == "trt"]) # 1
round(table(data2$death_flag21[data2$group == "trt"])/560 * 100, 2) # 0.18
## table sup - control unmatched ===============================================
### number individuals 
length(unique(ctr$eDWID)) # [1] 10294
### number of FEV1 measurements per individual
quantile2(table(ctr$eDWID)) # [1] 20.00 (13.00–29.00)
### total follow-up time (years) per individuals
head_rows <- tapply(seq_along(ctr$eDWID), ctr$eDWID, head, 1)
tail_rows <- tapply(seq_along(ctr$eDWID), ctr$eDWID, tail, 1)
quantile2(ctr$encounterage[tail_rows] - ctr$encounterage[head_rows]) # 4.18 (2.93–5.12)
remove(head_rows, tail_rows)
### year of the matching encounter
quantile2(ctr0$ReviewYear) # 2011.00 (2010.00–2012.00)
### sex
table2(ctr0$Gender) # 4701 (45.67%) 5593 (54.33%)
### age at the matching encounter
quantile2(ctr0$encounterage) # 15.25 (8.88–23.45)
### age cohort at the matching encounter
ctr0$age_cat <- cut(ctr0$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                    labels = c("<12", "[12,18)", "≥18"))
table2(ctr0$age_cat)
### annualized aztreonam
table2(ctr0$ann_aztreonam) # 8963 (87.07%) 1331 (12.93%)
### CF-related diabetes
table2(ctr0$cfrd_status2) # 8413 (81.73%) 385 (3.74%) 1496 (14.53%)
### annualized colistin
table2(ctr0$ann_colistin) # 9513 (92.41%) 781 (7.59%)
### annualized corticosteroids (oral)
table2(ctr0$ann_corticosteroids1) # 9631 (93.56%) 663 (6.44%)
### annualized corticosteroids (inhaled)
table2(ctr0$ann_corticosteroids2) # 7259 (70.52%) 3035 (29.48%)
### diagnosis age
quantile2(ctr0$Age_Diag) # 0.30 (0.00–1.40)
### annualized dornase alfa
table2(ctr0$ann_dornasealfa) # 1539 (14.95%) 8755 (85.05%)
### annualized aspergillus fumigatus
table2(ctr0$ann_fungalyeast1) # 8841 (85.88%) 1453 (14.12%)
### annualized ppFEV1
quantile2(ctr0$ann_GLI_FEV1_pct_predicted) # 83.92 (64.94–97.77)
### annualized ppFVC
quantile2(ctr0$ann_GLI_FVC_pct_predicted) # 92.80 (78.90–103.92)
### annualized ppFEV1/FVC
quantile2(ctr0$ann_GLI_FEV1FVC_pct_predicted) # 90.34 (80.66–97.09)
### annualized antifungals
table2(ctr0$ann_leukotriene_others3) # 9946 (96.62%) 348 (3.38%)
### annualized nontuberculous mycobacterium
table2(ctr0$ann_mycobacterials_nontub) # 9842 (95.61%) 452 (4.39%)
### annualized pseudomonas aeruginosa
table2(ctr0$ann_pseudomonasaeruginosa) # 5263 (51.13%) 5031 (48.87%)
### annualized methicillin-resistant staphylococcus aureus
table2(ctr0$ann_staphyl_type1) # 7597 (73.80%) 2697 (26.20%)
### annualizedi tobramycin (inhaled)
table2(ctr0$ann_tobi) # 5134 (49.87%) 5160 (50.13%)
## table sup - control matched =================================================
### number individuals 
length(unique(ctr2$eDWID)) # [1] 2800
### number of FEV1 measurements per individual
quantile2(table(ctr2$eDWID)) # 15.00 (9.00–23.00)
### total follow-up time (years) per individuals
head_rows <- tapply(seq_along(ctr2$eDWID), ctr2$eDWID, head, 1)
tail_rows <- tapply(seq_along(ctr2$eDWID), ctr2$eDWID, tail, 1)
quantile2(ctr2$encounterage[tail_rows] - ctr2$encounterage[head_rows]) # 3.18 (2.23–4.14)
remove(head_rows, tail_rows)
### year of the matching encounter
quantile2(ctr20$ReviewYear) # 2012.00 (2012.00–2013.00)
### sex
table2(ctr20$Gender) # 1253 (44.75%) 1547 (55.25%)
### age at the matching encounter
quantile2(ctr20$encounterage) # 15.10 (7.01–24.78)
### age cohort at the matching encounter
ctr20$age_cat <- cut(ctr20$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                     labels = c("<12", "[12,18)", "≥18"))
table2(ctr20$age_cat) # 1166 (41.64%) 466 (16.64%) 1168 (41.71%)
### annualized aztreonam
table2(ctr20$ann_aztreonam) # 2259 (80.68%) 541 (19.32%)
### CF-related diabetes
table2(ctr20$cfrd_status2) # 2287 (81.68%) 134 (4.79%) 379 (13.54%)
### annualized colistin
table2(ctr20$ann_colistin) # 2657 (94.89%) 143 (5.11%)
### annualized corticosteroids (oral)
table2(ctr20$ann_corticosteroids1) # 2645 (94.46%) 155 (5.54%)
### annualized corticosteroids (inhaled)
table2(ctr20$ann_corticosteroids2) # 1927 (68.82%) 873 (31.18%)
### diagnosis age
quantile2(ctr20$Age_Diag) # 0.40 (0.00–2.30)
### annualized dornase alfa
table2(ctr20$ann_dornasealfa) # 415 (14.82%) 2385 (85.18%)
### annualized aspergillus fumigatus
table2(ctr20$ann_fungalyeast1) # 2464 (88.00%) 336 (12.00%)
### annualized ppFEV1
quantile2(ctr20$ann_GLI_FEV1_pct_predicted) # 87.28 (68.91–100.75)
### annualized ppFVC
quantile2(ctr20$ann_GLI_FVC_pct_predicted) # 94.87 (82.31–106.24)
### annualized ppFEV1/FVC
quantile2(ctr20$ann_GLI_FEV1FVC_pct_predicted) # 91.70 (82.35–98.07)
### annualized antifungals
table2(ctr20$ann_leukotriene_others3) # 2683 (95.82%) 117 (4.18%)
### annualized nontuberculous mycobacterium
table2(ctr20$ann_mycobacterials_nontub) # 2640 (94.29%) 160 (5.71%)
### annualized pseudomonas aeruginosa
table2(ctr20$ann_pseudomonasaeruginosa) # 1509 (53.89%) 1291 (46.11%)
### annualized methicillin-resistant staphylococcus aureus
table2(ctr20$ann_staphyl_type1) # 2073 (74.04%) 727 (25.96%)
### annualizedi tobramycin (inhaled)
table2(ctr20$ann_tobi) # 1640 (58.57%) 1160 (41.43%)
## table sup - treated =========================================================
### number individuals 
length(unique(trt$eDWID)) # [1] 560
### number of FEV1 measurements per individual
quantile2(table(trt$eDWID)) # 15.00 (10.00–20.00)
### total follow-up time (years) per individuals
head_rows <- tapply(seq_along(trt$eDWID), trt$eDWID, head, 1)
tail_rows <- tapply(seq_along(trt$eDWID), trt$eDWID, tail, 1)
quantile2(trt$encounterage[tail_rows] - trt$encounterage[head_rows]) # 4.04 (3.20–4.47)
remove(head_rows, tail_rows)
### year of the matching encounter
quantile2(trt0$ReviewYear) # 2012.00 (2012.00–2012.00)
### sex
table2(trt0$Gender) # 251 (44.82%) 309 (55.18%)
### age at the matching encounter
quantile2(trt0$encounterage) # 15.62 (9.06–23.98)
### age cohort at the matching encounter
trt0$age_cat <- cut(trt0$age_0, breaks = c(-Inf, 12, 18, +Inf), right = FALSE, 
                    labels = c("<12", "[12,18)", "≥18"))
table2(trt0$age_cat) # 199 (35.54%) 133 (23.75%) 228 (40.71%)
### annualized aztreonam
table2(trt0$ann_aztreonam) # 450 (80.36%) 110 (19.64%)
### CF-related diabetes
table2(trt0$cfrd_status2) # 460 (82.14%) 26 (4.64%) 74 (13.21%)
### annualized colistin
table2(trt0$ann_colistin) # 530 (94.64%) 30 (5.36%)
### annualized corticosteroids (oral)
table2(trt0$ann_corticosteroids1) # 528 (94.29%) 32 (5.71%)
### annualized corticosteroids (inhaled)
table2(trt0$ann_corticosteroids2) # 388 (69.29%) 172 (30.71%)
### diagnosis age
quantile2(trt0$Age_Diag) # 0.60 (0.10–3.00)
### annualized dornase alfa
table2(trt0$ann_dornasealfa) # 85 (15.18%) 475 (84.82%)
### annualized aspergillus fumigatus
table2(trt0$ann_fungalyeast1) # 491 (87.68%) 69 (12.32%)
### annualized ppFEV1
quantile2(trt0$ann_GLI_FEV1_pct_predicted) # 88.44 (70.60–99.99)
### annualized ppFVC
quantile2(trt0$ann_GLI_FVC_pct_predicted) # 95.72 (82.51–106.24)
### annualized ppFEV1/FVC
quantile2(trt0$ann_GLI_FEV1FVC_pct_predicted) # 92.13 (82.46–98.12)
### annualized antifungals
table2(trt0$ann_leukotriene_others3) # 533 (95.18%) 27 (4.82%)
### annualized nontuberculous mycobacterium
table2(trt0$ann_mycobacterials_nontub) # 528 (94.29%) 32 (5.71%)
### annualized pseudomonas aeruginosa
table2(trt0$ann_pseudomonasaeruginosa) # 292 (52.14%) 268 (47.86%)
### annualized methicillin-resistant staphylococcus aureus
table2(trt0$ann_staphyl_type1) # 416 (74.29%) 144 (25.71%)
### annualizedi tobramycin (inhaled)
table2(trt0$ann_tobi) # 323 (57.68%) 237 (42.32%)
# between-subject analysis =====================================================
library(nlme)
library(splines)
m2a <- lme(GLI_FEV1_pct_predicted ~ (time * group) * age_cat,
           random = ~ 1 | subclass/eDWID,
           correlation = corExp(form =~ time | subclass/eDWID, nugget = TRUE),
           data = data2, method = "ML")
m2b <- lme(GLI_FEV1_pct_predicted ~ (ns(time, df = 2, B = c(0, 6)) * group) * age_cat,
           random = ~ 1 | subclass/eDWID,
           correlation = corExp(form =~ time | subclass/eDWID, nugget = TRUE),
           data = data2, method = "ML")
m2c <- lme(GLI_FEV1_pct_predicted ~ (ns(time, df = 3, B = c(0, 6)) * group) * age_cat,
           random = ~ 1 | subclass/eDWID,
           correlation = corExp(form =~ time | subclass/eDWID, nugget = TRUE),
           data = data2, method = "ML")
## table - model selection =====================================================
AIC(m2a, m2b, m2c)
BIC(m2a, m2b, m2c)
anova(m2a, m2b) # 5e-04
anova(m2a, m2c) # 0.015
## plot - model selection ======================================================
col <- rgb(0, 0, 0, alpha = 0.1)
get_est <- function(object, res = 100) {
  toy <- expand.grid(time = seq(0, 6, length.out = res),
                     age_cat = factor(c("<12", "[12,18)", "≥18"), 
                                      levels = c("<12", "[12,18)", "≥18")),
                     group = factor(c("ctr", "trt"), levels = c("ctr", "trt")))
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
est2a <- get_est(m2a)
est2b <- get_est(m2b)
est2c <- get_est(m2c)
png("manuscript/figures/model_selection2.png", width = 6 * 300, height = 8 * 300, res = 300)
{
  ylim <- range(est2a$lower, est2a$upper, est2b$lower, est2b$upper, est2c$lower, 
                est2c$upper)
  par(mfrow = c(3, 2), bg = "white", mar = c(5.1, 4.1, 4.1 - 2, 2.1))
  for(i in c("<12", "[12,18)", "≥18")) {
    for (j in c("ctr", "trt")) {
      bol <- est2a$age_cat == i & est2a$group == j
      plot(NA, xlim = c(0, 6), ylim = ylim, type = 'n', bty = "n", 
           xaxt = "n", yaxt = "n", 
           ylab = expression(ppFEV[1]), 
           xlab = "Time since matching encounter (years)")
      title <- paste0(ifelse(j == "ctr", "Control", "Treated"), ", ", i)
      mtext(title, side = 3, line = 0.5, 
            cex = 0.8)
      box(lwd = 0.5)
      axis(1, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5)
      axis(2, col = NA, col.ticks = 1, col.axis = 1, lwd.ticks = 0.5, las = 2)
      polygon(x = c(est2a$time[bol], rev(est2a$time[bol])), 
              y = c(est2a$lower[bol], rev(est2a$upper[bol])), col = col, border = NA)
      polygon(x = c(est2b$time[bol], rev(est2b$time[bol])), 
              y = c(est2b$lower[bol], rev(est2b$upper[bol])), col = col, border = NA)
      polygon(x = c(est2c$time[bol], rev(est2c$time[bol])), 
              y = c(est2c$lower[bol], rev(est2c$upper[bol])), col = col, border = NA)
      lines(est2a$time[bol], est2a$mean[bol])
      lines(est2b$time[bol], est2b$mean[bol])
      lines(est2c$time[bol], est2c$mean[bol])
      bol <- data2$age_cat == i & data2$group == j
      segments(x0 = data2$time[bol], x1 = data2$time[bol], y0 = par("usr")[3], 
               y1 = par("usr")[3] + 0.05 * diff(par("usr")[3:4]), col = col)
    }
  }
}
dev.off()
m2a <- update(m2a, method = "REML")
m2 <- m2a
m2 <- update(m2, fixed = GLI_FEV1_pct_predicted ~ -1 + time:group:age_cat + age_cat:group)
round(intervals(m2, which = "fixed")[[1]][, c(2, 1, 3)], 2) # estimates control vs treatment
m2 <- update(m2, fixed = GLI_FEV1_pct_predicted ~ -1 + (time * group) : age_cat + age_cat)
round(intervals(m2, which = "fixed")[[1]][, c(2, 1, 3)], 2) # ivacaftor effect
## effects plot ================================================================
png("manuscript/figures/effects2.png", width = 8 * 300, height = 5 * 300, res = 300)
{
  res <- 3
  toy <- expand.grid(time = seq(0, 6, length.out = res),
                     age_cat = factor(c("<12", "[12,18)", "≥18"), 
                                      levels = c("<12", "[12,18)", "≥18")),
                     group = factor(c("ctr", "trt"), levels = c("ctr", "trt")))
  frml <- formula(m2)[c(1, 3)]
  X <- model.matrix(frml, data = toy)
  betas <- fixef(m2)
  D <- vcov(m2)
  se <- apply(X, 1, function(x) sqrt(t(x) %*% D %*% x))
  toy$mean <- X %*% betas # mean
  toy$upper <- toy$mean + qnorm(0.975) * se # upper
  toy$lower <- toy$mean - qnorm(0.975) * se # lower
  ylim <- range(toy$lower, toy$upper)
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
    bol2 <- toy$age_cat == i & toy$group == "trt"
    polygon(x = c(toy$time[bol1], rev(toy$time[bol1])), 
            y = c(toy$lower[bol1], rev(toy$upper[bol1])), col = col, border = NA)
    polygon(x = c(toy$time[bol2], rev(toy$time[bol2])), 
            y = c(toy$lower[bol2], rev(toy$upper[bol2])), col = col, border = NA)
    lines(toy$time[bol1], toy$mean[bol1], lty = 2)
    lines(toy$time[bol2], toy$mean[bol2])
    legend("bottomleft", legend = c("Control", "On ivacaftor"),
           lwd = 1, lty = 2:1, bty = "n")
  }
}
dev.off()
## residuals ===================================================================
col <- rgb(0, 0, 0, alpha = 0.1)
res <- residuals(m2, type = "pearson")
fit <- fitted(m2)
png("manuscript/figures/residuals2.png", width = 10 * 300, height = 6 * 300, res = 300)
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
