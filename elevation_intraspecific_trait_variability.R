library(tidyverse)
library(vegan)
library(corrplot)
library(cati)
library(ggplot2)
library(cowplot)
library(patchwork)


####### dataset upload ########
plotxenv <- read.csv("/plotxenv.csv", row.names = 1) %>% 
  mutate(forest_type = as.factor(forest_type))

spxplot_inter <- read.csv("/spxplot_inter.csv", row.names = 1)

cwm_inter <- read.csv("/CWM_inter.csv", row.names = 1)
cwm_intra <- read.csv("/CWM_intra.csv", row.names = 1)

SES_FD_inter <- read.csv("/SES_FD_inter.csv", row.names = 1)
SES_FD_intra <- read.csv("/SES_FD_intra.csv", row.names = 1)

cwm_ITV <- cbind(
  cwm_intra %>% mutate(cwm_ITV_hveg = cwm_sp_hveg - cwm_inter$cwm_fx_hveg) %>% select(cwm_ITV_hveg),
  cwm_intra %>% mutate(cwm_ITV_LA = cwm_sp_la - cwm_inter$cwm_fx_la) %>% select(cwm_ITV_LA),
  cwm_intra %>% mutate(cwm_ITV_SLA = cwm_sp_sla - cwm_inter$cwm_fx_sla) %>% select(cwm_ITV_SLA),
  cwm_intra %>% mutate(cwm_ITV_LDMC = cwm_sp_ldmc - cwm_inter$cwm_fx_ldmc) %>% select(cwm_ITV_LDMC)
) #calculation of intraspecific trait variability for each trait

SES_FD_ITV <- cbind(
  SES_FD_intra %>% mutate(Rao_Hveg_ITV_SES = Rao_Hveg_sp_SES - SES_FD_inter$Rao_Hveg_fx_SES) %>% select(Rao_Hveg_ITV_SES),
  SES_FD_intra %>% mutate(Rao_LA_ITV_SES = Rao_LA_sp_SES - SES_FD_inter$Rao_LA_fx_SES) %>% select(Rao_LA_ITV_SES),
  SES_FD_intra %>% mutate(Rao_SLA_ITV_SES = Rao_SLA_sp_SES - SES_FD_inter$Rao_SLA_fx_SES) %>% select(Rao_SLA_ITV_SES),
  SES_FD_intra %>% mutate(Rao_LDMC_ITV_SES = Rao_LDMC_sp_SES - SES_FD_inter$Rao_LDMC_fx_SES) %>% select(Rao_LDMC_ITV_SES)
) #calculation of traspecific trait variability for each trait 

##################### correlation environmental variables ########################
####1. Standardization of environmental variables####
env_standardized <- plotxenv %>% 
  select(-c( belt, latitude, longitude, forest_type)) %>% 
  decostand(., method = "standardize") %>% 
  cbind(., plotxenv$forest_type) %>% 
  rename(forest_type = "plotxenv$forest_type") %>% 
  mutate(forest_type = as.factor(forest_type))

####2. Correlation and plot####
corr_db <- env_standardized %>% 
  rename("Elevation m" = altitude, "Slope °" = slope ,"Aspect 180°" = aspect_180,
         "Rocks %" = rocks_perc, "Stones %" = stones_perc,
         "Litter %" = litter_perc,
         "Tree %" = tree_perc, "Shrub %" = shrub_perc, "Herb %" = herb_perc,
         "Moss %" = moss_perc, "Deadwood %" = deadwood_perc, 
         "Soil depth cm" = soil_depth_med, "Canopy closure %" = canopy_closure)

file_path= "Correlation matrix.tiff" # Initialize file path

tiff(file=file_path,height = 250, width = 250, units = "mm", res = 300, compression = "lzw")

corrplot(cor(corr_db[,-c(14)], method = "pearson"), method = "circle", type = "full", order = "hclust", number.cex = 0.6, tl.col = "black" )

dev.off()     

################################### traitflexanova ######################################
#Lepš, J., de Bello, F., Šmilauer, P. and Doležal, J., 2011. Community trait response to environment: disentangling species turnover vs intraspecific trait variability effects. Ecography, 34(5), pp.856-863.
#####1. traitflex.anova ~ CWM ####

tfa_cwm_hveg <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, cwm_intra$cwm_sp_hveg, cwm_inter$cwm_fx_hveg, data=env_standardized)
tfa_cwm_sla <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, cwm_intra$cwm_sp_sla, cwm_inter$cwm_fx_sla, data=env_standardized)
tfa_cwm_la <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med  + forest_type, cwm_intra$cwm_sp_la, cwm_inter$cwm_fx_la, data=env_standardized)
tfa_cwm_ldmc <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, cwm_intra$cwm_sp_ldmc, cwm_inter$cwm_fx_ldmc, data=env_standardized)


#####2. traitflex.anova ~ SES_FD ####

tfa_ses_hveg <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, SES_FD_intra$Rao_Hveg_sp_SES, SES_FD_inter$Rao_Hveg_fx_SES, data=env_standardized)
tfa_ses_sla <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, SES_FD_intra$Rao_SLA_sp_SES, SES_FD_inter$Rao_SLA_fx_SES, data=env_standardized)
tfa_ses_la <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, SES_FD_intra$Rao_LA_sp_SES, SES_FD_inter$Rao_LA_fx_SES, data=env_standardized)
tfa_ses_ldmc <- traitflex.anova(~ altitude  + canopy_closure + soil_depth_med + forest_type, SES_FD_intra$Rao_LDMC_sp_SES, SES_FD_inter$Rao_LDMC_fx_SES, data=env_standardized)

#####3. ggplot traitflex.anova bargraph####

tabtr_cwm_hveg <- data.frame(t(tfa_cwm_hveg$RelSumSq*100))
tabtr_cwm_hveg_2 <- data.frame(tabtr_cwm_hveg$altitude......[c(1,2)], rownames(tabtr_cwm_hveg)[c(1,2)], rep("Plant Height", 2), rep(tabtr_cwm_hveg$altitude......[4], 2))
colnames(tabtr_cwm_hveg_2) <- c("val", "cat", "tr", "lim")

tabtr_cwm_sla <- data.frame(t(tfa_cwm_sla$RelSumSq*100))
tabtr_cwm_sla_2 <- data.frame(tabtr_cwm_sla$altitude......[c(1,2)], rownames(tabtr_cwm_sla)[c(1,2)], rep("SLA", 2), rep(tabtr_cwm_sla$altitude......[4], 2))
colnames(tabtr_cwm_sla_2) <- c("val", "cat", "tr", "lim")

tabtr_cwm_la <- data.frame(t(tfa_cwm_la$RelSumSq*100))
tabtr_cwm_la_2 <- data.frame(tabtr_cwm_la$altitude......[c(1,2)], rownames(tabtr_cwm_la)[c(1,2)], rep("Leaf Area", 2), rep(tabtr_cwm_la$altitude......[4], 2))
colnames(tabtr_cwm_la_2) <- c("val", "cat", "tr", "lim")

tabtr_cwm_ldmc <- data.frame(t(tfa_cwm_ldmc$RelSumSq*100))
tabtr_cwm_ldmc_2 <- data.frame(tabtr_cwm_ldmc$altitude......[c(1,2)], rownames(tabtr_cwm_ldmc)[c(1,2)], rep("LDMC", 2), rep(tabtr_cwm_ldmc$altitude......[4], 2))
colnames(tabtr_cwm_ldmc_2) <- c("val", "cat", "tr", "lim")

df_tfa_cwm <- rbind(tabtr_cwm_hveg_2, tabtr_cwm_la_2, tabtr_cwm_sla_2, tabtr_cwm_ldmc_2)

df_tfa_cwm$tr <- factor(df_tfa_cwm$tr, ordered = T, levels = c("LDMC", "SLA", "Leaf Area", "Plant Height"))

prop_exp_cwm <- ggplot(data = df_tfa_cwm, aes(y = tr, x = val, fill = cat)) +
  geom_bar(stat = "identity", color="darkgrey") +
  geom_errorbar(aes(x=lim,xmin=lim,xmax=lim))+
  labs(x = "Rel.sumsq.% explained", y = "") + 
  scale_fill_manual(values = c("white","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(color="black", size = 0.2)) +
  ggtitle("CWM") +
  theme(plot.title = element_text(hjust = 0.5, size = 25), legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 22), axis.title.x = element_text(size=24),
        legend.key.size = unit(2, 'cm'), legend.text = element_text(size=22))

tabtr_ses_hveg <- data.frame(t(tfa_ses_hveg$RelSumSq*100))
tabtr_ses_hveg_2 <- data.frame(tabtr_ses_hveg$altitude......[c(1,2)], rownames(tabtr_ses_hveg)[c(1,2)], rep("Plant Height", 2), rep(tabtr_ses_hveg$altitude......[4], 2))
colnames(tabtr_ses_hveg_2) <- c("val", "cat", "tr", "lim")

tabtr_ses_sla <- data.frame(t(tfa_ses_sla$RelSumSq*100))
tabtr_ses_sla_2 <- data.frame(tabtr_ses_sla$altitude......[c(1,2)], rownames(tabtr_ses_sla)[c(1,2)], rep("SLA", 2), rep(tabtr_ses_sla$altitude......[4], 2))
colnames(tabtr_ses_sla_2) <- c("val", "cat", "tr", "lim")

tabtr_ses_la <- data.frame(t(tfa_ses_la$RelSumSq*100))
tabtr_ses_la_2 <- data.frame(tabtr_ses_la$altitude......[c(1,2)], rownames(tabtr_ses_la)[c(1,2)], rep("Leaf Area", 2), rep(tabtr_ses_la$altitude......[4], 2))
colnames(tabtr_ses_la_2) <- c("val", "cat", "tr", "lim")

tabtr_ses_ldmc <- data.frame(t(tfa_ses_ldmc$RelSumSq*100))
tabtr_ses_ldmc_2 <- data.frame(tabtr_ses_ldmc$altitude......[c(1,2)], rownames(tabtr_ses_ldmc)[c(1,2)], rep("LDMC", 2), rep(tabtr_ses_ldmc$altitude......[4], 2))
colnames(tabtr_ses_ldmc_2) <- c("val", "cat", "tr", "lim")

df_tfa_ses <- rbind(tabtr_ses_hveg_2, tabtr_ses_sla_2, tabtr_ses_la_2, tabtr_ses_ldmc_2)

df_tfa_ses$tr <- factor(df_tfa_ses$tr, ordered = T, levels = c("LDMC", "SLA", "Leaf Area", "Plant Height"))

prop_exp_ses <- ggplot(data = df_tfa_ses, aes(y = tr, x = val, fill = cat)) +
  geom_bar(stat = "identity", color="darkgrey") +
  geom_errorbar(aes(x=lim,xmin=lim,xmax=lim))+
  labs(x = "Rel.sumsq.% explained", y = "") +
  guides(fill=guide_legend(title="")) + 
  scale_fill_manual(values = c("white","grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(color="black", size = 0.2)) +
  ggtitle("SES-FD") +
  theme(plot.title = element_text(hjust = 0.5, size = 25), axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 22), axis.title.x = element_text(size=24),
        legend.key.size = unit(2, 'cm'), legend.text = element_text(size=22))


variation_total <- plot_grid(prop_exp_cwm, prop_exp_ses, labels = c('A', 'B'),label_size = 25, ncol = 2)
ggsave(filename = "variation_total.tiff", plot = variation_total , width = 600, height = 300, units = "mm", dpi = 300, compression = "lzw")


############################### linear models ##############################
#####1.CWM#####

MyData_alt <- data.frame(
  altitude = seq(-1.6311, 1.7837, length = 100),
  soil_depth_med = 0,
  canopy_closure = 0,
  forest_type = levels(env_standardized$forest_type)[1])

##Hveg

mod_cwm_fx_hveg <- lm(cwm_inter$cwm_fx_hveg ~ altitude  + canopy_closure + soil_depth_med +  forest_type, data = env_standardized)

summary(mod_cwm_fx_hveg)

Pl_hveg_cwm_fx <- predict(mod_fx_hveg, newdata = MyData_alt, type = "response", se = TRUE)
pl_hveg_cwm_fx <- as.data.frame(do.call(cbind, lapply(Pl_hveg_cwm_fx, as.data.frame)))
colnames(pl_hveg_cwm_fx)[1] <- "fit_alt_cwm_hveg_fx"
colnames(pl_hveg_cwm_fx)[2] <- "se.fit_alt_cwm_hveg_fx"
colnames(pl_hveg_cwm_fx)[3] <- "residual.scale_alt_cwm_hveg_fx"
colnames(pl_hveg_cwm_fx)[4] <- "forest_type_fx"

MyData_alt_hveg_fx <- cbind(MyData_alt, pl_hveg_cwm_fx)

altitude_ns <- MyData_alt_hveg$altitude * sd(plotxenv$altitude) + mean(plotxenv$altitude)

MyData_alt_hveg_fx<- cbind(MyData_alt_hveg_fx, altitude_ns)

levels(plotxenv$forest_type)[levels(plotxenv$forest_type) == "beech"] <- "Beech"
levels(plotxenv$forest_type)[levels(plotxenv$forest_type) == "mixed"] <- "Mixed"

ggplot_cwm_turn_hveg <- ggplot(data = MyData_alt_hveg) + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_inter$cwm_fx_hveg, shape = forest_type),  size=2.5)+
  geom_smooth(aes(x=altitude_ns, y = fit_alt_cwm_hveg_fx), color="black", se=F, size = 0.5) +                              
  geom_smooth(aes(x=altitude_ns, y= (fit_alt_cwm_hveg_fx +1.96 * se.fit_alt_cwm_hveg_fx )), color="black", se=F,linetype=2)+
  geom_smooth(aes(x=altitude_ns, y= (fit_alt_cwm_hveg_fx -1.96 * se.fit_alt_cwm_hveg_fx )), color="black", se=F,linetype=2)+
  scale_shape_manual(values=c(2,3))+
  labs(x = "", y = "CWM Plant Height", shape = "Forest Type")+
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Turnover")

mod_cwm_sp_hveg <- lm(cwm_intra$cwm_sp_hveg ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_sp_hveg)

Pl_hveg_cwm_sp <- predict(mod_sp_hveg, newdata = MyData_alt, type = "response", se = TRUE)
pl_hveg_cwm_sp <- as.data.frame(do.call(cbind, lapply(Pl_hveg_cwm_sp, as.data.frame)))
colnames(pl_hveg_cwm_sp)[1] <- "fit_alt_cwm_hveg_sp"
colnames(pl_hveg_cwm_sp)[2] <- "se.fit_alt_cwm_hveg_sp"
colnames(pl_hveg_cwm_sp)[3] <- "residual.scale_alt_cwm_hveg_sp"
colnames(pl_hveg_cwm_sp)[4] <- "forest_type_sp"

MyData_alt_hveg_sp <- cbind(MyData_alt, pl_hveg_cwm_sp, altitude_ns)

ggplot_cwm_tot_hveg <- ggplot(data = MyData_alt_hveg_sp) + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_intra$cwm_sp_hveg, shape = forest_type),  size=2.5)+
  geom_smooth(aes(x=altitude_ns, y = fit_alt_cwm_hveg_sp), color="black", se=F, size = 0.5) +                              
  geom_smooth(aes(x=altitude_ns, y= (fit_alt_cwm_hveg_sp +1.96 * se.fit_alt_cwm_hveg_sp )), color="black", se=F,linetype=2)+
  geom_smooth(aes(x=altitude_ns, y= (fit_alt_cwm_hveg_sp -1.96 * se.fit_alt_cwm_hveg_sp )), color="black", se=F,linetype=2)+
  scale_shape_manual(values=c(2,3))+
  labs(x = "", y = "CWM Plant Height", shape = "Forest Type")+
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Total")

mod_cwm_itv_hveg <- lm(cwm_ITV$cwm_ITV_hveg ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_itv_hveg)

ggplot_cwm_itv_hveg <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_ITV$cwm_ITV_hveg, shape=forest_type), size=2.5)+
  labs(x = "", y = "CWM Plant Height", shape = "Forest Type")+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("ITV")

##CWM_LA
mod_cwm_fx_la <- lm(cwm_inter$cwm_fx_la ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_fx_la)

ggplot_cwm_turn_la <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_inter$cwm_fx_la, shape=forest_type), size=2.5)+
  labs(x = "", y = "CWM Leaf area", shape = "Forest Type")+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Turnover")

mod_cwm_sp_la <- lm(cwm_intra$cwm_sp_la ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_sp_la)

ggplot_cwm_tot_la <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_intra$cwm_sp_la, shape=forest_type), size=2.5)+
  labs(x = "", y = "CWM Leaf area", shape = "Forest Type")+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Total")

mod_cwm_itv_la <- lm(cwm_ITV$cwm_ITV_LA ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_itv_la)

ggplot_cwm_itv_la <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_ITV$cwm_ITV_LA, shape=forest_type), size=2.5)+
  labs(x = "", y = "CWM Leaf area", shape = "Forest Type")+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("ITV")

##CWM_SLA
mod_cwm_fx_sla <- lm(cwm_inter$cwm_fx_sla ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_fx_sla)

ggplot_cwm_turn_sla <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_inter$cwm_fx_sla, shape=forest_type), size=2.5)+
  labs(x = "", y = "CWM SLA", shape = "Forest Type") + 
  scale_shape_manual(values = c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_cwm_sp_sla <- lm(cwm_intra$cwm_sp_sla ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_sp_sla)

ggplot_cwm_tot_sla <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_intra$cwm_sp_sla, shape=forest_type), size=2.5)+
  labs(x = "", y = "CWM SLA", shape = "Forest Type")+ 
  scale_shape_manual(values = c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_cwm_itv_sla <- lm(cwm_ITV$cwm_ITV_SLA ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_itv_sla)

ggplot_cwm_itv_sla <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_ITV$cwm_ITV_SLA, shape=forest_type), size=2.5) +
  labs(x = "Elevation", y = "CWM SLA", shape = "Forest Type") + 
  scale_shape_manual(values = c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

##CWM_LDMC
mod_cwm_fx_ldmc <- lm(cwm_inter$cwm_fx_ldmc ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_fx_ldmc)

Pl_ldmc_cwm_fx <- predict(mod_cwm_fx_ldmc, newdata = MyData_alt, type = "response", se = TRUE)
pl_ldmc_cwm_fx <- as.data.frame(do.call(cbind, lapply(Pl_ldmc_cwm_fx, as.data.frame)))
colnames(pl_ldmc_cwm_fx)[1] <- "fit_alt_cwm_ldmc_fx"
colnames(pl_ldmc_cwm_fx)[2] <- "se.fit_alt_cwm_ldmc_fx"
colnames(pl_ldmc_cwm_fx)[3] <- "residual.scale_alt_cwm_ldmc_fx"
colnames(pl_ldmc_cwm_fx)[4] <- "forest_type_fx"

MyData_alt_ldmc_fx <- cbind(MyData_alt, pl_ldmc_cwm_fx, altitude_ns)

ggplot_cwm_turn_ldmc <- ggplot(data = MyData_alt_ldmc_fx) + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_inter$cwm_fx_ldmc, shape = forest_type), size=2.5)+
  geom_smooth(aes(x=altitude_ns, y = fit_alt_cwm_ldmc_fx), color="black", se=F) +                              
  geom_smooth(aes(x=altitude_ns, y= (fit_alt_cwm_ldmc_fx +1.96 * se.fit_alt_cwm_ldmc_fx )), color="black", se=F,linetype=2)+
  geom_smooth(aes(x=altitude_ns, y= (fit_alt_cwm_ldmc_fx -1.96 * se.fit_alt_cwm_ldmc_fx )), color="black", se=F,linetype=2)+
  labs(x = "", y = "CWM LDMC", shape = "Forest Type")+
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_cwm_sp_ldmc <- lm(cwm_intra$cwm_sp_ldmc ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_sp_ldmc)

ggplot_cwm_tot_ldmc <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_intra$cwm_sp_ldmc, shape = forest_type), size=2.5)+
  labs(x = "", y = "CWM LDMC", shape = "Forest Type")+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_cwm_itv_ldmc <- lm(cwm_ITV$cwm_ITV_LDMC ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_cwm_itv_ldmc)

ggplot_cwm_itv_ldmc <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=cwm_ITV$cwm_ITV_LDMC, shape = forest_type), size=2.5)+
  labs(x = "Elevation (m)", y = "CWM LDMC", shape = "Forest Type")+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

#plot constructions: Hveg + LDMC; LA + SLA

ggplot_cwm_hveg_ldmc <- (ggplot_cwm_turn_hveg | ggplot_cwm_itv_hveg | ggplot_cwm_tot_hveg)/(ggplot_cwm_turn_ldmc | ggplot_cwm_itv_ldmc | ggplot_cwm_tot_ldmc)
ggplot_CWM_H_LDMC <- ggplot_cwm_hveg_ldmc + plot_layout(guides = "collect")
ggsave(filename = "ggplot_cwm_hveg_ldmc.tiff", plot = ggplot_CWM_H_LDMC, width = 400, height = 340, units = "mm", dpi = 300, compression = "lzw")

ggplot_cwm_la_sla <- (ggplot_cwm_turn_la | ggplot_cwm_itv_la | ggplot_cwm_tot_la)/(ggplot_cwm_turn_sla | ggplot_cwm_itv_sla | ggplot_cwm_tot_sla)
ggplot_CWM_LA_SLA <- ggplot_cwm_la_sla + plot_layout(guides = "collect")
ggsave(filename = "ggplot_cwm_la_sla.tiff", plot = ggplot_CWM_LA_SLA, width = 400, height = 340, units = "mm", dpi = 300, compression = "lzw")


#####2.SES-FD####

##Hveg
mod_ses_fd_fx_hveg <- lm(SES_FD_inter$Rao_Hveg_fx_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_fx_hveg)

ggplot_ses_fd_turn_hveg <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_inter$Rao_Hveg_fx_SES, shape = forest_type),  size=2.5)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ylim(-3, 3)+
  labs(x = "", y = "SES-FD Plant Height", shape = "Forest Type")+
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Turnover")

mod_ses_fd_sp_hveg <- lm(SES_FD_intra$Rao_Hveg_sp_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_sp_hveg)

ggplot_ses_fd_tot_hveg <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_intra$Rao_Hveg_sp_SES, shape = forest_type), size=2.5)+
  labs(x = "", y = "SES-FD Plant Height", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  geom_smooth(data = plotxenv, aes(x = altitude, y = SES_FD_intra$Rao_Hveg_sp_SES), method = "lm", se = FALSE, color = "blue") +
  ylim(-3, 3)+
  scale_shape_manual(values=c(2,3)) + 
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Total") 

mod_ses_fd_itv_hveg <- lm(SES_FD_ITV$Rao_Hveg_ITV_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_itv_hveg)

ggplot_ses_fd_itv_hveg <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_ITV$Rao_Hveg_ITV_SES, shape = forest_type), size=2.5)+
  geom_smooth(data = plotxenv, aes(x = altitude, y = SES_FD_ITV$Rao_Hveg_ITV_SES), method = "lm", se = FALSE, color = "blue") +
  labs(x = "", y = "SES-FD Plant Height", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3)+ 
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("ITV")

##LA
mod_ses_fd_fx_la <- lm(SES_FD_inter$Rao_LA_fx_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_fx_la)

ggplot_ses_fd_turn_la <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_inter$Rao_LA_fx_SES, shape = forest_type), size=2.5)+
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x = "", y = "SES-FD Leaf Area", shape = "Forest Type") +
  ylim(-3, 3) +
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Turnover")

mod_ses_fd_sp_la <- lm(SES_FD_intra$Rao_LA_sp_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_sp_la)

ggplot_ses_fd_tot_la <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_intra$Rao_LA_sp_SES, shape = forest_type), size=2.5)+
  labs(x = "", y = "SES-FD Leaf Area", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3) +
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("Total")

mod_ses_fd_itv_la <- lm(SES_FD_ITV$Rao_LA_ITV_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_itv_la)

ggplot_ses_fd_itv_la <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_ITV$Rao_LA_ITV_SES, shape = forest_type), size=2.5)+
  labs(x = "", y = "SES-FD Leaf Area", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3) +
  scale_shape_manual(values=c(2,3)) +
  theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
         axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
         panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("ITV")

##SLA
mod_ses_fd_fx_sla <- lm(SES_FD_inter$Rao_SLA_fx_SES ~ altitude + soil_depth_med + canopy_closure + forest_type, env_standardized)

summary(mod_ses_fd_fx_sla)

ggplot_ses_fd_turn_sla <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_inter$Rao_SLA_fx_SES, shape = forest_type), size=2.5)+
  geom_hline(yintercept=0, linetype="dashed")+
  labs(x = "", y = "SES-FD SLA", shape = "Forest Type")+
  ylim(-3, 3) +
  scale_shape_manual(values=c(2,3)) +
  theme( panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
         axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
         panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_ses_fd_sp_sla <- lm(SES_FD_intra$Rao_SLA_sp_SES ~ altitude + soil_depth_med + canopy_closure + forest_type, env_standardized)

summary(mod_ses_fd_sp_sla)

ggplot_ses_fd_tot_sla <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_intra$Rao_SLA_sp_SES, shape = forest_type), size=2.5)+
  labs(x = "", y = "SES-FD SLA", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3) +
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_ses_fd_itv_sla <- lm(SES_FD_ITV$Rao_SLA_ITV_SES ~ altitude + soil_depth_med + canopy_closure + forest_type, env_standardized)

summary(mod_ses_fd_itv_sla)

ggplot_ses_fd_itv_sla <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_ITV$Rao_SLA_ITV_SES, shape = forest_type), size=2.5)+
  geom_smooth(data = plotxenv, aes(x = altitude, y = SES_FD_ITV$Rao_SLA_ITV_SES), method = "lm", se = FALSE, color = "blue") +
  labs(x = "Elevation", y = "SES-FD SLA", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3) +
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")


##LDMC
mod_ses_fd_fx_ldmc <- lm(SES_FD_inter$Rao_LDMC_fx_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_fx_ldmc)

ggplot_ses_fd_turn_ldmc <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_inter$Rao_LDMC_fx_SES, shape = forest_type), size=2.5)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_smooth(data = plotxenv, aes(x = altitude, y = SES_FD_inter$Rao_LDMC_fx_SES), method = "lm", se = FALSE, color = "blue") +
  ylim(-3, 3)+
  scale_shape_manual(values=c(2,3)) +
  labs(x = "", y = "SES-FD LDMC", shape = "Forest Type")+
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_ses_fd_sp_ldmc <- lm(SES_FD_intra$Rao_LDMC_sp_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_sp_ldmc)

ggplot_ses_fd_tot_ldmc <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_intra$Rao_LDMC_sp_SES, shape = forest_type), size=2.5)+
  geom_smooth(data = plotxenv, aes(x = altitude, y = SES_FD_intra$Rao_LDMC_sp_SES), method = "lm", se = FALSE, color = "blue") +
  labs(x = "", y = "SES-FD LDMC", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3)+
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=18),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

mod_ses_fd_itv_ldmc <- lm(SES_FD_ITV$Rao_LDMC_ITV_SES ~ altitude  + canopy_closure + soil_depth_med + forest_type, data = env_standardized)

summary(mod_ses_fd_itv_ldmc)

ggplot_ses_fd_itv_ldmc <- ggplot() + 
  geom_point(data = plotxenv, aes(x= altitude, y=SES_FD_ITV$Rao_LDMC_ITV_SES, shape = forest_type), size=2.5)+
  geom_smooth(data = plotxenv, aes(x = altitude, y = SES_FD_ITV$Rao_LDMC_ITV_SES), method = "lm", se = FALSE, color = "blue") +
  labs(x = "Elevation (m)", y = "SES-FD LDMC", shape = "Forest Type")+ 
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-3, 3)+
  scale_shape_manual(values=c(2,3)) +
  theme(panel.grid.major = element_blank(), axis.title.y=element_text(size=18),
        axis.text.y=element_text(colour = "black", size=16),axis.text.x=element_text(colour = "black", size=16),
        panel.grid.minor = element_blank(),axis.title.x=element_text(size=26),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=22)) +
  theme(legend.key.size = unit(2, 'cm'), legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  ggtitle ("")

#plot constructions: Hveg + LDMC; LA + SLA

ggplot_ses_fd_hveg_ldmc <- (ggplot_ses_fd_turn_hveg | ggplot_ses_fd_itv_hveg | ggplot_ses_fd_tot_hveg)/(ggplot_ses_fd_turn_ldmc | ggplot_ses_fd_itv_ldmc | ggplot_ses_fd_tot_ldmc)
ggplot_SES_FD_H_LDMC <- ggplot_ses_fd_hveg_ldmc + plot_layout(guides = "collect")
ggsave(filename = "ggplot_ses_fd_hveg_ldmc.tiff", plot = ggplot_SES_FD_H_LDMC, width = 400, height = 340, units = "mm", dpi = 300, compression = "lzw")

ggplot_ses_fd_la_sla <- (ggplot_ses_fd_turn_la | ggplot_ses_fd_itv_la | ggplot_ses_fd_tot_la)/(ggplot_ses_fd_turn_sla | ggplot_ses_fd_itv_sla | ggplot_ses_fd_tot_sla)
ggplot_SES_FD_LA_SLA <- ggplot_ses_fd_la_sla + plot_layout(guides = "collect")
ggsave(filename = "ggplot_ses_fd_la_sla.tiff", plot = ggplot_SES_FD_LA_SLA, width = 400, height = 340, units = "mm", dpi = 300, compression = "lzw")


###################### t.test - convergence/random/divergence assembly determination##############################

Hveg_fx <- t.test(SES_FD_inter$Rao_Hveg_fx_SES, mu = 0)

LA_fx <- t.test(SES_FD_inter$Rao_LA_fx_SES, mu = 0)

SLA_fx <- t.test(SES_FD_inter$Rao_SLA_fx_SES, mu = 0)

LDMC_fx <- t.test(SES_FD_inter$Rao_LDMC_fx_SES, mu = 0)

Hveg_sp <- t.test(SES_FD_intra$Rao_Hveg_sp_SES, mu = 0)

LA_sp <- t.test(SES_FD_intra$Rao_LA_sp_SES, mu = 0)

SLA_sp <- t.test(SES_FD_intra$Rao_SLA_sp_SES, mu = 0)

LDMC_sp <- t.test(SES_FD_intra$Rao_LDMC_sp_SES, mu = 0)

Hveg_itv <- t.test(SES_FD_ITV$Rao_Hveg_ITV_SES, mu = 0)

LA_itv <- t.test(SES_FD_ITV$Rao_LA_ITV_SES, mu = 0)

SLA_itv <- t.test(SES_FD_ITV$Rao_SLA_ITV_SES, mu = 0)

LDMC_itv <- t.test(SES_FD_ITV$Rao_LDMC_ITV_SES, mu = 0)


########################## species overlap ##########################
db_sp_inter <- read.csv("/spxplot_inter.csv", row.names = 1)

DCA <- decorana(log1p(t(db_sp_inter)))
DCA

plotxenv$elevation <- plotxenv$altitude

ef <- envfit(DCA ~ elevation, data = plotxenv, perm = 999)
tiff(file="dca_plot.tiff",height = 250, width = 250, units = "mm", res = 300, compression = "lzw")
ordiplot(DCA, display = "sp", type = "text", col = "black")
plot(ef, col = "red")

dev.off()

