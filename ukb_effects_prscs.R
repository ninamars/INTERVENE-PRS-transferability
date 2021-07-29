
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# PRS-CS associations in UK Biobank

# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# South Asian and African / Caribbean

library(data.table)
library(dplyr)
library(stringi)
library(pROC)
library(ROCR)
library(ggplot2)

d <- fread("phenos_intervene_20201026.txt")

cad_prs <- fread("prscs_cad.add.160614.website_UKB_nonwhite_PRS_20210615.txt")
names(cad_prs) <- c("eid", "cadprs")
d <- merge(d, cad_prs, by = "eid")

brcanc_prs <- fread("prscs_michailidoubreastcancerall_UKB_nonwhite_PRS_20210615.txt")
names(brcanc_prs) <- c("eid", "brcancprs")
d <- merge(d, brcanc_prs, by = "eid")

prcanc_prs <- fread("prscs_meta_v3_onco_euro_overall_ChrAll_1_release_UKB_nonwhite_PRS_20210615.txt")
names(prcanc_prs) <- c("eid", "prcancprs")
d <- merge(d, prcanc_prs, by = "eid")

t2d_prs <- fread("prscs_METAANALYSIS_DIAGRAM_SE1_UKB_nonwhite_PRS_20210615.txt")
names(t2d_prs) <- c("eid", "t2dprs")
d <- merge(d, t2d_prs, by = "eid")

# -------------------------------------------------------------------------

# Coronary artery disease

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "chd", "eth2", "cadprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$chd)
summary(d$chd_age[d$eid %in% temp$eid])  # Age at diagnosis exists only for cases, applies for all diseases
sd(d$chd_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("chd ~ scale(cadprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

#exp(coefficients(fit)[2])
#exp(coefficients(fit)[2] - 1.96*summary(fit)$coefficients[2,2])
#exp(coefficients(fit)[2] + 1.96*summary(fit)$coefficients[2,2])

# -------------------------------------------------------------------------

# Type 2 diabetes

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
        select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "diabetes", "eth2", "t2dprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$diabetes)
summary(d$diabetes_age[d$eid %in% temp$eid])
sd(d$diabetes_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("diabetes ~ scale(t2dprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

# -------------------------------------------------------------------------

# Breast cancer

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
      select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "brcanc", "eth2", "brcancprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$brcanc)
summary(d$brcanc_age[d$eid %in% temp$eid])
sd(d$brcanc_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("brcanc ~ scale(brcancprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

# -------------------------------------------------------------------------

# Prostate cancer

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "prcanc", "eth2", "prcancprs"); dim(temp)
temp <- na.omit(temp);dim(temp) 
table(temp$prcanc)
summary(d$prcanc_age[d$eid %in% temp$eid])
sd(d$prcanc_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("prcanc ~ scale(prcancprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# EUR ancestry


d <- fread("phenos_intervene_20201026.txt")

cad_prs <- fread("cad.add.160614.website_UKB_PRS_20201105.txt")
names(cad_prs) <- c("eid", "cadprs")
d <- merge(d, cad_prs, by = "eid")

brcanc_prs <- fread("prscs_michailidoubreastcancerall_UKB_PRS_20201028.txt")
names(brcanc_prs) <- c("eid", "brcancprs")
d <- merge(d, brcanc_prs, by = "eid")

prcanc_prs <- fread("onco_euro_overall_ChrAll_1_release_UKB_PRS_20201105.txt")
names(prcanc_prs) <- c("eid", "prcancprs")
d <- merge(d, prcanc_prs, by = "eid")

t2d_prs <- fread("METAANALYSIS_DIAGRAM_SE1_UKB_PRS_20201105.txt")
names(t2d_prs) <- c("eid", "t2dprs")
d <- merge(d, t2d_prs, by = "eid")

# -------------------------------------------------------------------------

# Coronary artery disease

covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d %>% select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "chd", "eth2", "cadprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$chd)
summary(d$chd_age[d$eid %in% temp$eid])
sd(d$chd_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("chd ~ scale(cadprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

#exp(coefficients(fit)[2])
#exp(coefficients(fit)[2] - 1.96*summary(fit)$coefficients[2,2])
#exp(coefficients(fit)[2] + 1.96*summary(fit)$coefficients[2,2])

# -------------------------------------------------------------------------

# Type 2 diabetes

covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d %>% select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "diabetes", "eth2", "t2dprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$diabetes)
summary(d$diabetes_age[d$eid %in% temp$eid])
sd(d$diabetes_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("diabetes ~ scale(t2dprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

# -------------------------------------------------------------------------

# Breast cancer

covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d %>% select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "brcanc", "eth2", "brcancprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$brcanc)
summary(d$brcanc_age[d$eid %in% temp$eid])
sd(d$brcanc_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("brcanc ~ scale(brcancprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

# -------------------------------------------------------------------------

# Prostate cancer

covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d %>% select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "prcanc", "eth2", "prcancprs"); dim(temp)
temp <- na.omit(temp);dim(temp)
table(temp$prcanc)
summary(d$prcanc_age[d$eid %in% temp$eid])
sd(d$prcanc_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("prcanc ~ scale(prcancprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

exp(coefficients(fit)[2])
exp(coefficients(fit)[2] - 1.96*summary(fit)$coefficients[2,2])
exp(coefficients(fit)[2] + 1.96*summary(fit)$coefficients[2,2])

# -------------------------------------------------------------------------
