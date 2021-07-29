
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

# Plot UK Biobank effects across ancestries for PRS-CS vs LDpred vs small scores

# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(readxl)

d <- read_excel("Figure2_results_detailed.xlsx")
d$endpoint <- factor(d$endpoint, levels = c("CAD", "T2D", "Breast cancer", "Prostate cancer"))
d$prs_type <- factor(d$prs_type, levels = c("Limited-variant PRS", "LDpred", "PRS-CS"))

plot_effects <- function(temp = temp,
                         ancestry = ancestry) {
  
  return(ggplot(temp,
                aes(x = endpoint, y = or, group = prs_type, fill = prs_type))+
           geom_hline(yintercept = 1.0, lwd = 0.2, col = "gray4")+
           geom_point(size = 2.9, position = position_dodge(width=.55), col = "gray15", pch = 22)+
           geom_errorbar(aes(ymin=lower, ymax=upper), width=0, lwd = 0.35, position=position_dodge(.55), col = "gray10")+
           theme_bw()+ theme(panel.grid.major.x = element_blank())+
           coord_cartesian(ylim = c(0.9, 3))+
           labs(y = "OR per SD (95% CI)", x = "", title = ancestry)+
           theme(plot.title=element_text(size=10,face="bold"))+
           theme(axis.title=element_text(size=9))+
           theme(plot.margin = margin(0.1,0,0,0.5, "cm"))+
           theme(text = element_text(size = 10))+
           theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))+
           theme(legend.title = element_blank(), legend.text=element_text(size=8.5))+
           theme(legend.position = "right")+
           scale_fill_manual(values = c("#6D6D68", "#ECD623", "#FF7C4B"))
  )
}

p_eur <- plot_effects(temp = d %>% filter(group == "European"),
                      ancestry = "UK Biobank, European")
p_eur

p_sa <- plot_effects(temp = d %>% filter(group == "South Asian"),
                     ancestry = "UK Biobank, South Asian")
p_sa

p_afr <- plot_effects(temp = d %>% filter(group == "African / Caribbean"),
                      ancestry = "UK Biobank, African / Caribbean")
p_afr

ggarrange(p_eur, p_sa, p_afr, nrow = 1, ncol = 3,
          common.legend = T, legend = "right")
