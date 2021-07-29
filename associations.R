
######## ANALYSES PERFORMED SIMILARLY ACROSS BIOBANKS AND COHORTS

# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# PRS-disease associations in UK Biobank
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# South Asian and African / Caribbean

library(data.table)
library(dplyr)
library(stringi)
library(pROC)
library(ROCR)
library(ggplot2)

d <- fread("phenos_intervene_20201026.txt")

#### Replace the input files with relevant PRS file, here showing PRS-CS PRSs

cad_prs <- fread("prscs_cad.add.160614.website_UKB_nonwhite_PRS.txt")
names(cad_prs) <- c("eid", "cadprs")
d <- merge(d, cad_prs, by = "eid")

brcanc_prs <- fread("prscs_michailidoubreastcancerall_UKB_nonwhite_PRS.txt")
names(brcanc_prs) <- c("eid", "brcancprs")
d <- merge(d, brcanc_prs, by = "eid")

prcanc_prs <- fread("prscs_meta_v3_onco_euro_overall_ChrAll_1_release_UKB_nonwhite_PRS.txt")
names(prcanc_prs) <- c("eid", "prcancprs")
d <- merge(d, prcanc_prs, by = "eid")

t2d_prs <- fread("prscs_METAANALYSIS_DIAGRAM_SE1_UKB_nonwhite_PRS.txt")
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

cad_prs <- fread("cad.add.160614.website_UKB_PRS.txt")
names(cad_prs) <- c("eid", "cadprs")
d <- merge(d, cad_prs, by = "eid")

brcanc_prs <- fread("prscs_michailidoubreastcancerall_UKB_PRS.txt")
names(brcanc_prs) <- c("eid", "brcancprs")
d <- merge(d, brcanc_prs, by = "eid")

prcanc_prs <- fread("onco_euro_overall_ChrAll_1_release_UKB_PRS.txt")
names(prcanc_prs) <- c("eid", "prcancprs")
d <- merge(d, prcanc_prs, by = "eid")

t2d_prs <- fread("METAANALYSIS_DIAGRAM_SE1_UKB_PRS.txt")
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




# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# PRS-disease associations in FinnGen across settlement regions
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plyr)

data <- fread(...)

data <- data[data$AGE_AT_DEATH_OR_NOW >= 18,]

late <- c("North Karelia", "North Ostrobothnia", "Pohjois-Savo", "Kainuu") 
mixed <- c("Kymenlaakso", "Lapland", "South Karelia", "Etelä-Savo", "Central Finland", "Etelä-Karjala")
early <- c("Kanta-Häme", "Uusimaa", "South Ostrobothnia", "Varsinais-Suomi", "Pirkanmaa",
          "Päijät-Häme", "Ostrobothnia", "Satakunta", "Central Ostrobothnia", "Etelä-Pohjanmaa")

data$settlement <- with(data, ifelse(regionofbirthname %in% early, "Early settlement",
                                ifelse(regionofbirthname %in% mixed, "Borderline",
                                  ifelse(regionofbirthname %in% late, "Late settlement", NA))))
table(data$settlement, useNA = "always")
table(data$regionofbirthname, data$settlement, useNA = "always")
table(data$regionofbirthname, useNA = "always")

set <- c("Early settlement", "Late settlement", "Borderline")


# Define logistic regression formula with scaling PRS to mean zero and unit variance by settelement region,
# adjusting for age, sex, batches and 10PCs
formula <- formula(...)

res <- matrix(NA, nrow = length(set), ncol = 12)

for(ii in 1:3) {
  
  alue <- set[ii]
  res[ii,1] <- alue
  data_alue <- data[data$settlement == alue,]

  temp <- addmargins(table(data_alue$ENDPOINT)) # replace with relevant endpoint
  res[ii,2] <- temp[[1]]
  res[ii,3] <- temp[[2]]
  res[ii,4] <- temp[[3]]
  
  fit <- glm(formula, data = data_alue, family = "binomial")
  summary(fit)
  
  ... # Add relevant information to results file and save file
  
}









