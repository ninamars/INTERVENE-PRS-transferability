
library(data.table)
library(dplyr)
library(base)
library(readxl)
library(ggplot2)
library(GGally)
library(survival)
library(tidyr)
library(ggrepel)
library(meta)

#   -----------------------------------------------------------------------

# Meta-analyse effects across EUR ancestry cohorts:

met <- read_excel("Figure1_results_detailed.xlsx") %>%
  filter(ancestry == "European")
met

d <- c()

for(endpoint in unique(met$endpoint)) {
    print(paste("Starting:", endpoint))
    
    res <- metagen(log(met$or[met$endpoint == endpoint]),
                   met$se[met$endpoint == endpoint], sm = "OR")
    d <- rbind(d, c(endpoint, exp(c(res$TE.fixed, res$lower.fixed, res$upper.fixed))))
    
}

d <- data.frame(d)
names(d) <- c("endpoint", "or", "lower", "upper")
d$col2 <- rep("turquoise4", 4)
d$order <- rep(5, 4)
d$group <- "European"

write.table(d, "Figure1_results_eur_meta_fixed.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(d)

#   -----------------------------------------------------------------------

######### Panel A

d <- read_excel("Figure1_results_detailed.xlsx") %>%
  filter(ancestry != "European") %>% select(c("endpoint", "or", "lower", "upper", "col2", "order", "group"))
d <- rbind(d, fread("Figure1_results_eur_meta_fixed.txt"))

d$endpoint <- factor(d$endpoint, levels = c("CAD", "T2D", "Breast cancer", "Prostate cancer"))

d <- d %>% 
  ungroup() %>% 
  arrange(order) %>% 
  group_by(endpoint) %>% 
  mutate(.r = row_number())
d

p1 <- ggplot(d %>% group_by(col2), aes(x = reorder(group, .r), y = or))+
  geom_hline(yintercept = 1.0, lwd = 0.2, col = "gray4")+
  geom_point(size = 3, pch = 15, color = d$col2, alpha = 0.9) +
  theme_bw()+ theme(panel.grid.major.x = element_blank())+
  facet_grid(. ~ endpoint, scales = "free")+
  coord_cartesian(ylim = c(0.5, 2.25))+
  labs(y = "OR (95% CI)", x = "")+
  theme(plot.margin = margin(0.5,0.5,0 ,2, "cm"))+
  geom_errorbar(aes(ymin=(lower), ymax=(upper)), width=0, lwd = 0.35,
                position=position_dodge(.9))+
  theme(text = element_text(size = 10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
  theme(aspect.ratio = 1.2)
  
p1

pdf("Fig1_panelA.pdf")
plot(p1)
dev.off()
rm(d)

######### Panel B

d <- read_excel("Figure1_results_detailed.xlsx") %>%
  filter(ancestry == "European") %>% select(c("endpoint", "or", "lower", "upper", "col2", "order", "group"))
d$endpoint <- factor(d$endpoint, levels = c("CAD", "T2D", "Breast cancer", "Prostate cancer"))

d <- d %>% 
  ungroup() %>% 
  arrange(or) %>% 
  mutate(.r = row_number())
d

p2 <- ggplot(d %>% group_by(col2), aes(x = reorder(group, .r), y = or))+
  geom_hline(yintercept = 1.0, lwd = 0.2, col = "gray4")+
  geom_point(size = 3, pch = 15, color = d$col2, alpha = 0.9) +
  theme_bw()+ theme(panel.grid.major.x = element_blank())+
  facet_grid(. ~ endpoint, scales = "free")+
  coord_cartesian(ylim = c(0.5, 2.25))+
  labs(y = "OR (95% CI)", x = "")+
  theme(plot.margin = margin(0.5,0.5,0 ,2, "cm"))+
  geom_errorbar(aes(ymin=(lower), ymax=(upper)), width=0, lwd = 0.35,
                position=position_dodge(.9))+
  theme(text = element_text(size = 10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
  theme(aspect.ratio = 1.2)
p2

pdf("Fig1_panelB.pdf")
plot(p2)
dev.off()
rm(d)


######### Panel C

d <- read_excel("Figure1_data_panelC.xlsx")

d$endpoint <- factor(d$endpoint, levels = c("CAD", "T2D", "Breast cancer", "Prostate cancer"))
d$settlement <- factor(d$settlement, levels = c("Early settlement", "Borderline", "Late settlement"))

p3 <- ggplot(d, aes(x = settlement, y = or))+
  geom_hline(yintercept = 1.0, lwd = 0.2, col = "gray4")+
  geom_point(size = 3, pch = 15, color = d$col) +
  theme_bw()+ theme(panel.grid.major.x = element_blank())+
  facet_grid(. ~ endpoint, scales = "free")+
  coord_cartesian(ylim = c(0.5, 2.25))+
  labs(y = "OR (95% CI)", x = "")+
  theme(plot.margin = margin(0.5,0.5,0 ,2, "cm"))+
  geom_errorbar(aes(ymin=(lower), ymax=(upper)), width=0, lwd = 0.35,
                position=position_dodge(.9))+
  theme(text = element_text(size = 10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
  theme(aspect.ratio = 1.2)

p3

pdf("Fig1_panelC.pdf")
plot(p3)
dev.off()
