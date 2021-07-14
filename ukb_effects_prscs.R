
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# PRS-CS associations in UK Biobank

# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# South Asian and African / Caribbean

/fs/projects/finngen/incoming_data/nmars/apps/R-3.4.2/bin/R

library(data.table)
library(dplyr)
library(stringi)
library(pROC)
library(ROCR)
library(ggplot2)

d <- fread("/fs/projects/ukbb/nmars/phenos/phenos_intervene_20201026.txt")

cad_prs <- fread("/fs/projects/ukbb/nmars/R6_prscs_cad.add.160614.website_UKB_nonwhite_PRS_20210615.txt")
names(cad_prs) <- c("eid", "cadprs")
d <- merge(d, cad_prs, by = "eid")

brcanc_prs <- fread("/fs/projects/ukbb/nmars/finngen_R5_michailidoubreastcancerall_UKB_nonwhite_PRS_20210615.txt")
names(brcanc_prs) <- c("eid", "brcancprs")
d <- merge(d, brcanc_prs, by = "eid")

prcanc_prs <- fread("/fs/projects/ukbb/nmars/R6_prscs_meta_v3_onco_euro_overall_ChrAll_1_release_UKB_nonwhite_PRS_20210615.txt")
names(prcanc_prs) <- c("eid", "prcancprs")
d <- merge(d, prcanc_prs, by = "eid")

t2d_prs <- fread("/fs/projects/ukbb/nmars/R6_prscs_finngen_R6_METAANALYSIS_DIAGRAM_SE1_UKB_nonwhite_PRS_20210615.txt")
names(t2d_prs) <- c("eid", "t2dprs")
d <- merge(d, t2d_prs, by = "eid")

table(d$eth2)
tapply(d$baselineage, d$eth2, summary)
tapply(d$baselineage, d$eth2, sd, na.rm = T)
tapply(d$sex, d$eth2, function(x){ prop.table(table(x)) })

# -------------------------------------------------------------------------

# Coronary artery disease

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
                select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "chd", "eth2", "cadprs"); dim(temp)
temp <- na.omit(temp);dim(temp) # 7628 SA; 7618 Black/Caribbean
table(temp$chd)
summary(d$chd_age[d$eid %in% temp$eid])
sd(d$chd_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("chd ~ scale(cadprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

# -------------------------------------------------------------------------

# Type 2 diabetes

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "sex", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
        select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "diabetes", "eth2", "t2dprs"); dim(temp)
temp <- na.omit(temp);dim(temp) # 7260 SA; 7346 Black/Caribbean
table(temp$diabetes)
summary(d$diabetes_age[d$eid %in% temp$eid])
sd(d$diabetes_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("diabetes ~ scale(t2dprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)

#exp(coefficients(fit)[2])
#exp(coefficients(fit)[2] - 1.96*summary(fit)$coefficients[2,2])
#exp(coefficients(fit)[2] + 1.96*summary(fit)$coefficients[2,2])

# -------------------------------------------------------------------------

# Breast cancer

ancestry <- "Black / Caribbean"
ancestry <- "South Asian"

covs <- paste(c("baselineage", "batch", paste0("C", 1:10)), collapse = " + ")
temp <- d[d$eth2 == ancestry,] %>%
      select("eid", "baselineage", "sex", "batch", paste0("C", 1:10), "brcanc", "eth2", "brcancprs"); dim(temp)
temp <- na.omit(temp);dim(temp) # 3514 SA; 4342 Black/Caribbean
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
temp <- na.omit(temp);dim(temp) # 4d114 SA; 3276 Black/Caribbean
table(temp$prcanc)
summary(d$prcanc_age[d$eid %in% temp$eid])
sd(d$prcanc_age[d$eid %in% temp$eid], na.rm = T)

formula <- formula(paste0("prcanc ~ scale(prcancprs) + ", covs))
fit <- glm(formula, data = temp, family = "binomial")
summary(fit)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# EUR ancestry

/fs/projects/finngen/incoming_data/nmars/apps/R-3.4.2/bin/R

library(data.table)
library(dplyr)
library(stringi)
library(pROC)
library(ROCR)
library(ggplot2)

d <- fread("/fs/projects/ukbb/nmars/phenos/phenos_intervene_20201026.txt")

cad_prs <- fread("/fs/projects/ukbb/nmars/cad.add.160614.website_UKB_R6_PRS_20201105.txt")
names(cad_prs) <- c("eid", "cadprs")
d <- merge(d, cad_prs, by = "eid")

brcanc_prs <- fread("/fs/projects/ukbb/nmars/prscs_michailidoubreastcancerall_UKB_PRS_20201028.txt")
names(brcanc_prs) <- c("eid", "brcancprs")
d <- merge(d, brcanc_prs, by = "eid")

prcanc_prs <- fread("/fs/projects/ukbb/nmars/onco_euro_overall_ChrAll_1_release_UKB_R6_PRS_20201105.txt")
names(prcanc_prs) <- c("eid", "prcancprs")
d <- merge(d, prcanc_prs, by = "eid")

t2d_prs <- fread("/fs/projects/ukbb/nmars/METAANALYSIS_DIAGRAM_SE1_UKB_R6_PRS_20201105.txt")
names(t2d_prs) <- c("eid", "t2dprs")
d <- merge(d, t2d_prs, by = "eid")

dim(d) # 343676

table(d$eth2)
tapply(d$baselineage, d$eth2, summary)
tapply(d$baselineage, d$eth2, sd, na.rm = T)
tapply(d$sex, d$eth2, function(x){ prop.table(table(x)) })

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
