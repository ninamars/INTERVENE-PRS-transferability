
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# LDpred effects across all tested fractions of causal variants in UKB

# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

library(data.table)
library(dplyr)
library(stringi)
library(pROC)
library(ROCR)
library(ggpubr)
library(ggplot2)

thresholds <- c("p1.0000e-04", "p3.0000e-04", "p1.0000e-03", "p3.0000e-03", "p1.0000e-02", "p3.0000e-02", "p1.0000e-01", "p3.0000e-01", "p1.0000e+00", "inf")


# -------------------------------------------------------------------------

# CAD:

d <- fread("phenos_intervene_20201026.txt")

# EUR:
cad_prs <- fread("cad1000g_UKB_white_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)
# Non-EUR:
cad_prs <- fread("cad1000g_UKB_nonwhite_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)

d <- merge(d, cad_prs, by = "eid")
covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")

# EUR:
cad <- select(d, c("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "chd", "eth2", thresholds)); dim(cad)

# African / Caribbean:
cad <- d %>% filter(eth2 == "African / Caribbean") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "chd", "eth2", thresholds); dim(cad)
# South Asian
cad <- d %>% filter(eth2 == "South Asian") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "chd", "eth2", thresholds); dim(cad)

# cad <- na.omit(cad);dim(cad) # does not include NAs
names(cad)[17:26] <- c(paste0("prs", 1:10))
table(cad$chd)
summary(d$chd_age[d$eid %in% cad$eid])
sd(d$chd_age[d$eid %in% cad$eid], na.rm = T)

res <- matrix(NA, nrow = 10, ncol = 5)

for(ii in 1:10) {
  formula <- formula(paste0("chd ~ scale(prs", ii,") + ", covs))
  fit <- glm(formula, data = cad, family = "binomial")
  res[ii,] <- c(paste0("prs", ii), summary(fit)$coef[2,])
}

write.table(res, paste0("cad_sa_or_per_sd_by_p.txt"), # replace sa with bc for African / Caribbean ancestry
            quote = F, col.names = T,
            sep = "\t", na = "", row.names = F)

# -------------------------------------------------------------------------

# T2D:

d <- fread("phenos_intervene_20201026.txt")

# EUR:
t2d_prs <- fread("dm2with1000g_UKB_white_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)
# Non-EUR:
t2d_prs <- fread("dm2with1000g_UKB_nonwhite_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)

d <- merge(d, t2d_prs, by = "eid")
covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")

# EUR:
t2d <- select(d, c("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "diabetes", "eth2", thresholds)); dim(t2d)

# African / Caribbean:
t2d <- d %>% filter(eth2 == "African / Caribbean") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "diabetes", "eth2", thresholds); dim(t2d)
# South Asian
t2d <- d %>% filter(eth2 == "South Asian") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "diabetes", "eth2", thresholds); dim(t2d)

t2d <- na.omit(t2d);dim(t2d) # excludes T1D patients and those diagnosed before age 18
names(t2d)[17:26] <- c(paste0("prs", 1:10))
table(t2d$diabetes)
summary(d$diabetes_age[d$eid %in% t2d$eid])
sd(d$diabetes_age[d$eid %in% t2d$eid], na.rm = T)

res <- matrix(NA, nrow = 10, ncol = 5)

for(ii in 1:10) {
  formula <- formula(paste0("diabetes ~ scale(prs", ii,") + ", covs))
  fit <- glm(formula, data = t2d, family = "binomial")
  res[ii,] <- c(paste0("prs", ii), summary(fit)$coef[2,])
}

write.table(res, paste0("t2d_sa_or_per_sd_by_p.txt"), # replace sa with bc for African / Caribbean ancestry
            quote = F, col.names = T,
            sep = "\t", na = "", row.names = F)

# -------------------------------------------------------------------------

# Breast cancer:

d <- fread("phenos_intervene_20201026.txt")

# EUR:
brcanc_prs <- fread("michailidoubreastcancerall1000g_UKB_white_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)
# Non-EUR:
brcanc_prs <- fread("michailidoubreastcancerall1000g_UKB_nonwhite_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)

d <- merge(d, brcanc_prs, by = "eid")
covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ") # No sex as covariate

# EUR:
brcanc <- select(d, c("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "brcanc", "eth2", thresholds)); dim(brcanc)

# African / Caribbean:
brcanc <- d %>% filter(eth2 == "African / Caribbean") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "brcanc", "eth2", thresholds); dim(brcanc)
# South Asian
brcanc <- d %>% filter(eth2 == "South Asian") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "brcanc", "eth2", thresholds); dim(brcanc)

brcanc <- na.omit(brcanc);dim(brcanc) 
names(brcanc)[17:26] <- c(paste0("prs", 1:10))
table(brcanc$brcanc)
summary(d$brcanc_age[d$eid %in% brcanc$eid])
sd(d$brcanc_age[d$eid %in% brcanc$eid], na.rm = T)

res <- matrix(NA, nrow = 10, ncol = 5)

for(ii in 1:10) {
  formula <- formula(paste0("brcanc ~ scale(prs", ii,") + ", covs)) # men are NA in brcanc variable
  fit <- glm(formula, data = brcanc, family = "binomial")
  res[ii,] <- c(paste0("prs", ii), summary(fit)$coef[2,])
}

write.table(res, paste0("brcanc_sa_or_per_sd_by_p.txt"), # replace sa with bc for African / Caribbean ancestry
            quote = F, col.names = T,
            sep = "\t", na = "", row.names = F)

# -------------------------------------------------------------------------

# Prostate cancer:

d <- fread("phenos_intervene_20201026.txt")

# EUR:
prcanc_prs <- fread("schumacherprostatecancer1000g_UKB_white_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)
# Non-EUR:
prcanc_prs <- fread("schumacherprostatecancer1000g_UKB_nonwhite_LDpred_GRS_v2.txt") %>%
          select(c("IID", thresholds)) %>% rename(eid = IID)

d <- merge(d, prcanc_prs, by = "eid")
covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ") # No sex as covariate

# EUR:
prcanc <- select(d, c("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "prcanc", "eth2", thresholds)); dim(prcanc)

# African / Caribbean:
prcanc <- d %>% filter(eth2 == "African / Caribbean") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "prcanc", "eth2", thresholds); dim(prcanc)
# South Asian
prcanc <- d %>% filter(eth2 == "South Asian") %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "prcanc", "eth2", thresholds); dim(prcanc)

prcanc <- na.omit(prcanc);dim(prcanc) 
names(prcanc)[17:26] <- c(paste0("prs", 1:10))
table(prcanc$prcanc)
summary(d$prcanc_age[d$eid %in% prcanc$eid])
sd(d$prcanc_age[d$eid %in% prcanc$eid], na.rm = T)

res <- matrix(NA, nrow = 10, ncol = 5)

for(ii in 1:10) {
  formula <- formula(paste0("prcanc ~ scale(prs", ii,") + ", covs)) # women are NA in prcanc variable
  fit <- glm(formula, data = prcanc, family = "binomial")
  res[ii,] <- c(paste0("prs", ii), summary(fit)$coef[2,])
}

write.table(res, paste0("prcanc_sa_or_per_sd_by_p.txt"), # replace sa with bc for African / Caribbean ancestry
            quote = F, col.names = T,
            sep = "\t", na = "", row.names = F)

# -------------------------------------------------------------------------

# Plot result


diseases <- c("cad", "t2d", "brcanc", "prcanc")

thresholds <- c("p1.0000e-04", "p3.0000e-04", "p1.0000e-03",
                "p3.0000e-03", "p1.0000e-02", "p3.0000e-02",
                "p1.0000e-01", "p3.0000e-01", "p1.0000e+00", "inf")

all <- c()
for(disease in diseases) {
  files <- list.files(pattern = disease)
  all <- rbind(all, do.call(rbind, lapply(files, function(f) { cbind(fread(f), f) })))
}

# Assign labels for plots:
all$disease <-  factor(c(rep("Coronary artery disease", 30),
                         rep("Type 2 diabetes", 30),
                         rep("Breast cancer", 30),
                         rep("Prostate cancer", 30)))
all$disease_short <-  c(rep(diseases[1], 30),
                        rep(diseases[2], 30),
                        rep(diseases[3], 30),
                        rep(diseases[4], 30))
all$group <- factor(rep(c(rep("African / Caribbean", 10),
                          rep("European", 10),
                          rep("South Asian", 10)), 4),
                    levels = c("European", "South Asian", "African / Caribbean"))
all$threshold <- factor(rep(thresholds, 12), levels = thresholds)


plot_effects <- function(temp = temp,
                         disease = disease,
                         disease_short = disease_short,
                         face = face) {
  
  return(ggplot(temp, aes(x = threshold, y = exp(V2), group = group, fill = group))+
           # facet_grid(.~disease, scale = "free", space = "free", rows = 4)+
           geom_hline(yintercept = 1.0, lwd = 0.2, col = "gray4")+
           geom_point(size = 2.75, position = position_dodge(width=.55), col = "gray15", pch = 22)+
           geom_errorbar(aes(ymin=exp(V2-(1.96*V3)), ymax=exp(V2+(1.96*V3))), width=0, lwd = 0.25,
                         position=position_dodge(.55), col = "gray10")+
           theme_bw()+ theme(panel.grid.major.x = element_blank())+
           coord_cartesian(ylim = c(0.9, 2.5))+
           labs(y = "OR per SD (95% CI)", x = "", title = disease)+
           theme(plot.title=element_text(size=10,face="bold"))+
           theme(axis.title=element_text(size=9))+
           theme(plot.margin = margin(0.1,0,0,0.5, "cm"))+
           theme(text = element_text(size = 10))+
           theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))+
           theme(legend.title = element_blank(), legend.text=element_text(size=8.5))+
           theme(legend.position = "right")+
           scale_fill_manual(values = c("turquoise4", "slateblue4", "indianred"))+
           theme(axis.text.x = element_text(face = face)))
}

p_chd <- plot_effects(disease = "Coronary artery disease",
                      disease_short = "cad",
                      temp = all %>% filter(disease == "Coronary artery disease"),
                      face = c(rep("plain", 3),rep("bold", 1), rep("plain", 6)))

p_t2d <- plot_effects(disease = "Type 2 diabetes",
                      disease_short = "t2d",
                      temp = all %>% filter(disease == "Type 2 diabetes"),
                      face = c(rep("plain", 3),rep("bold", 1), rep("plain", 6)))

p_brcanc <- plot_effects(disease = "Breast cancer",
                         disease_short = "brcanc",
                         temp = all %>% filter(disease == "Breast cancer"),
                         face = c(rep("plain", 5),rep("bold", 1), rep("plain", 4)))

p_prcanc <- plot_effects(disease = "Prostate cancer",
                         disease_short = "prcanc",
                         temp = all %>% filter(disease == "Prostate cancer"),
                         face = c(rep("plain", 4),rep("bold", 1), rep("plain", 5)))

ggarrange(p_chd, p_t2d, p_brcanc, p_prcanc, nrow = 2, ncol = 2,
          common.legend = T, legend = "right")


