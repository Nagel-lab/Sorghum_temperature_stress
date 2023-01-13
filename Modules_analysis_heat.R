#######################
# Modules analysis heat
#######################

rm(list = ls())
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/heat_analysis/")

library(reshape2)
library(tidyr)
library(stringr)
library(doBy)
library(ggplot2)
library(dplyr)

# Plot the eigengenes for each module
#####################################
MEs_RTX_heat <- read.table("MEs_RTX_heat.txt", header = T)
MEs_Macia_heat <- read.table("MEs_Macia_heat.txt", header = T)
MEs_RTX_heat$sample <- row.names(MEs_RTX_heat)
MEs_Macia_heat$sample <- row.names(MEs_Macia_heat)

MEs_RTX_heat <- melt(MEs_RTX_heat, id.vars = "sample")
MEs_Macia_heat <- melt(MEs_Macia_heat, id.vars = "sample")
names(MEs_RTX_heat)[2] <- "module"
names(MEs_Macia_heat)[2] <- "module"
MEs_RTX_heat$Analysis <- "Modules_RTX"
MEs_Macia_heat$Analysis <- "Modules_Macia"

MEs_heat <- rbind(MEs_RTX_heat, MEs_Macia_heat)

# Separate the sample column
MEs_heat <- MEs_heat %>% separate(sample, into = c("Exp","Genotype", "Time", "Response"))
MEs_heat$Genotype <- str_replace_all(MEs_heat$Genotype, 'RTX', 'RTX430')
MEs_heat$Time <- str_replace_all(MEs_heat$Time, 'T', '')
MEs_heat$Response <- str_replace_all(MEs_heat$Response, 'C1', '30.1')
MEs_heat$Response <- str_replace_all(MEs_heat$Response, 'C2', '30.2')
MEs_heat$Response <- str_replace_all(MEs_heat$Response, 'C3', '30.3')
MEs_heat$Response <- str_replace_all(MEs_heat$Response, 'H1', '42.1')
MEs_heat$Response <- str_replace_all(MEs_heat$Response, 'H2', '42.2')
MEs_heat$Response <- str_replace_all(MEs_heat$Response, 'H3', '42.3')
MEs_heat <- MEs_heat %>% separate(Response, into = c("Temperature", "Replicate"))

# calculate mean and SD by module
#--------------------------------
MEs_heat_mean <- summary_by(.~Analysis+Genotype+Time+Temperature+module, keep.names = T,FUN = function(x){mean(x,na.rm=T)} , data= MEs_heat)
MEs_heat_sd <- summary_by(.~Analysis+Genotype+Time+Temperature+module, keep.names = T,FUN = function(x){sd(x,na.rm=T)} , data= MEs_heat)
names(MEs_heat_sd)[6] <- "SD"
MEs_heat_mean <- merge.data.frame(MEs_heat_mean, MEs_heat_sd, by = c("Analysis","Genotype","Time","Temperature","module"))
names(MEs_heat_mean)[6] <- "Mean"

# Plot the eigengene by module and by analysis
#---------------------------------------------
str(MEs_heat_mean)
MEs_heat_mean$Analysis <- as.factor(MEs_heat_mean$Analysis)
MEs_heat_mean$Genotype <- factor(MEs_heat_mean$Genotype, levels = c("RTX430","Macia"))
MEs_heat_mean$Time <- as.numeric(as.character(MEs_heat_mean$Time))
MEs_heat_mean$Temperature <- as.factor(MEs_heat_mean$Temperature)

# Replace the default modules names by numbers from 1 to X (e.g. RTX-H1 for module 1 in RTX430 in the heat experiment)
Analysis_RTX <- MEs_heat_mean[MEs_heat_mean$Analysis == "Modules_RTX",]
modules_names_RTX <- paste("RTX-H",c(seq(1,7,1),0,seq(8,17,1)), sep = "")
module <- c("MEsaddlebrown","MEantiquewhite4","MEturquoise", "MEbrown", "MEgreen", "MEtan",
            "MEdarkorange2", "MEgrey", "MEbrown4", "MEdarkred", "MEpaleturquoise", 
            "MElightpink4", "MEdarkturquoise", "MEmagenta", "MEdarkorange", "MEmediumpurple3", "MEbisque4", "MElightcyan")
names_modules_RTX <- as.data.frame(cbind(modules_names_RTX, module))
Analysis_RTX <- left_join(Analysis_RTX, names_modules_RTX, by = "module")
Analysis_RTX$modules_names_RTX <- factor(Analysis_RTX$modules_names_RTX, levels = modules_names_RTX)

setwd("C:/Users/tbonnot/Documents/sorghum/Paper/Plots/")

pdf(file = "MEs_profiles_RTX_heat_figure.pdf", width = 10, height = 7)
# the module grey is removed because it corresponds to genes that were unsigned to any modules
ggplot(Analysis_RTX[!Analysis_RTX$module == "MEgrey",], aes(x = Time, y = Mean, group= Temperature)) +
  geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.08, color = NA)+
  geom_point(data=Analysis_RTX[!Analysis_RTX$module == "MEgrey",],aes(shape = Temperature, colour = Temperature),
             position=position_dodge(width=0.1), size = 3)+
  geom_line(data=Analysis_RTX[!Analysis_RTX$module == "MEgrey",],position=position_dodge(width=0.1),
            aes(linetype = Temperature, colour = Temperature))+
  scale_shape_manual(values=c(16,17))+
  scale_linetype_manual(values=c("solid", "blank"))+
  geom_errorbar(data=Analysis_RTX[!Analysis_RTX$module == "MEgrey",],aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                position=position_dodge(width=0.1))+
  scale_color_manual(values=c("#666666","#FF6600"))+
  scale_x_continuous(breaks = seq(0,15,3))+
  theme_bw()+
  facet_wrap(~modules_names_RTX)+
  xlab("Time of day (ZT)")+
  ylab("Normalized transcript abundance of module eigengenes")+
  theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
        axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12))

dev.off()

Analysis_Macia <- MEs_heat_mean[MEs_heat_mean$Analysis == "Modules_Macia",]
modules_names_Macia <- paste("Macia-H",seq(1,20,1), sep = "")
module <- c("MElightpink4","MEthistle1","MEbrown","MEgrey60", "MEsaddlebrown", "MEdarkorange", "MEcoral1",
            "MEgreen", "MEmediumpurple3", "MEsienna3", "MEdarkolivegreen", "MEcyan", 
             "MEplum", "MEblack", "MEivory", "MEdarkslateblue", "MEhoneydew1", "MElightcyan1",
            "MEdarkmagenta","MEgreenyellow")
names_modules_Macia <- as.data.frame(cbind(modules_names_Macia, module))
Analysis_Macia <- left_join(Analysis_Macia, names_modules_Macia, by = "module")
Analysis_Macia$modules_names_Macia <- factor(Analysis_Macia$modules_names_Macia, levels = modules_names_Macia)

pdf(file = "MEs_profiles_Macia_heat_figure.pdf", width = 10, height = 7)

ggplot(Analysis_Macia, aes(x = Time, y = Mean, group= Temperature)) +
  geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.08, color = NA)+
  geom_point(data=Analysis_Macia,aes(shape = Temperature, colour = Temperature),
             position=position_dodge(width=0.1), size = 3)+
  geom_line(data=Analysis_Macia,position=position_dodge(width=0.1),
            aes(linetype = Temperature, colour = Temperature))+
  scale_shape_manual(values=c(16,17))+
  scale_linetype_manual(values=c("solid", "blank"))+
  geom_errorbar(data=Analysis_Macia,aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                position=position_dodge(width=0.1))+
  scale_color_manual(values=c("#666666","#FF6600"))+
  scale_x_continuous(breaks = seq(0,15,3))+
  theme_bw()+
  facet_wrap(~modules_names_Macia)+
  xlab("Time of day (ZT)")+
  ylab("Normalized transcript abundance of module eigengenes")+
  theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
        axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12))

dev.off()

# Plot a specific module
#-----------------------
ggplot(Analysis_Macia[Analysis_Macia$module == "MEblack",], aes(x = Time, y = Mean, group= Temperature)) +
  geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.08, color = NA)+
  geom_point(data=Analysis_Macia[Analysis_Macia$module == "MEblack",],
             aes(shape = Temperature, colour = Temperature),
             position=position_dodge(width=0.1), size = 3)+
  geom_line(data=Analysis_Macia[Analysis_Macia$module == "MEblack",],
            position=position_dodge(width=0.1),
            aes(linetype = Temperature, colour = Temperature))+
  scale_shape_manual(values=c(16,17))+
  scale_linetype_manual(values=c("solid", "blank"))+
  geom_errorbar(data=Analysis_Macia[Analysis_Macia$module == "MEblack",],
                aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                position=position_dodge(width=0.1))+
  scale_color_manual(values=c("#666666","#FF6600"))+
  scale_x_continuous(breaks = seq(0,15,3))+
  theme_bw()+
  #facet_wrap(~module)+
  xlab("Time of day (ZT)")+
  ylab("Normalized transcript abundance of module eigengenes")+
  theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
        axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12))


# Enrichment of DEGs in the different modules
#############################################
setwd("C:/Users/tbonnot/Documents/sorghum/LRT_new/")
load(file = "dds_LRT_Temperature_Genotype_heat.RData")
load(file = "dds_LRT_Genotype.heat.RData")

setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/heat_analysis/")
dds_LRT_Genotype.heat.signif <- dds_LRT_Genotype.heat[dds_LRT_Genotype.heat$padj < 0.05,]
dds_LRT_Genotype.heat.signif <- na.omit(dds_LRT_Genotype.heat.signif)

# module counts table Macia
mod_RTX_heat <- read.table("modules_RTX_heat.txt", header = T)
mod_Macia_heat <- read.table("modules_Macia_heat.txt", header = T)

Macia_module_count <- data.frame(table(mod_Macia_heat$moduleColors))
names(Macia_module_count) <- c("Module","Freq_all")
Macia_module_count$Non_all <- sum(Macia_module_count$Freq_all)-(Macia_module_count$Freq_all)

# Look at the 5989 genes with a genotype effect and see if they are enriched in any module in Macia
#--------------------------------------------------------------------------------------------------
DEGs_to_look <- dds_LRT_Genotype.heat.signif
DEGs_to_look$ID <- row.names(DEGs_to_look)
names(mod_Macia_heat)[1] <- "ID"
DEGs_to_look_module <- merge.data.frame(mod_Macia_heat, DEGs_to_look, by = "ID")
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
Macia_module_count$Proportion_all <- (Macia_module_count$Freq_all)/(sum(Macia_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(Macia_module_count, DEGs_to_look_module_count, by = "Module")

# Calculate fold enrichment
module_count_stat$Fold_enrichment <- (module_count_stat$Proportion_DEG)/(module_count_stat$Proportion_all)

# Calculate significance using Fisher's exact tests
list.y.var <- unique(module_count_stat$Module)

Sub = list()
result = list()

for(i in list.y.var)
{
  sub = module_count_stat[module_count_stat$Module == i, c(5,6,2,3)]
  names(sub) <- rep(c("N","non_N"),2)
  sub = rbind(sub[,1:2],sub[,3:4])
  sub = as.matrix(sub)
  fisher <- fisher.test(sub)
  result[[i]] <- fisher.test(sub)
  Sub[[i]] <- fisher$p.value
}

fisher <- do.call("rbind", lapply(Sub, as.data.frame))
fisher <- data.frame(row.names(fisher),fisher)
names(fisher) <- c("Module", "Fisher_pval")

enrich_geno <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_geno$group <- "Geno"


# Look at TFs and see if they are enriched in any module in Macia
#----------------------------------------------------------------
setwd("C:/Users/tbonnot/Documents/sorghum/PlantTFDB/")
TF_sorghum <- read.table("Sbi_TF_list.txt", header = T)
DEGs_to_look <- TF_sorghum$Gene_ID
DEGs_to_look_module <- filter(mod_Macia_heat, ID %in% DEGs_to_look)
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
Macia_module_count$Proportion_all <- (Macia_module_count$Freq_all)/(sum(Macia_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(Macia_module_count, DEGs_to_look_module_count, by = "Module")

# Calculate fold enrichment
module_count_stat$Fold_enrichment <- (module_count_stat$Proportion_DEG)/(module_count_stat$Proportion_all)

# Calculate significance using Fisher's exact tests
list.y.var <- unique(module_count_stat$Module)

Sub = list()
result = list()

for(i in list.y.var)
{
  sub = module_count_stat[module_count_stat$Module == i, c(5,6,2,3)]
  names(sub) <- rep(c("N","non_N"),2)
  sub = rbind(sub[,1:2],sub[,3:4])
  sub = as.matrix(sub)
  fisher <- fisher.test(sub)
  result[[i]] <- fisher.test(sub)
  Sub[[i]] <- fisher$p.value
}

fisher <- do.call("rbind", lapply(Sub, as.data.frame))
fisher <- data.frame(row.names(fisher),fisher)
names(fisher) <- c("Module", "Fisher_pval")

enrich_Macia_TF <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_Macia_TF$group <- "Macia_TF"


# Look at cycling genes and see if they are enriched in any module in Macia
#---------------------------------------------------------------------------
setwd("C:/Users/tbonnot/Documents/sorghum/Phase_enrichment/")
Cycling_sorghum <- read.csv2("Phase_reference_Lai_rounded.csv", header = T, sep = ",")
DEGs_to_look <- Cycling_sorghum$ï..AGI
DEGs_to_look_module <- filter(mod_Macia_heat, ID %in% DEGs_to_look)
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
Macia_module_count$Proportion_all <- (Macia_module_count$Freq_all)/(sum(Macia_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(Macia_module_count, DEGs_to_look_module_count, by = "Module")

# Calculate fold enrichment
module_count_stat$Fold_enrichment <- (module_count_stat$Proportion_DEG)/(module_count_stat$Proportion_all)

# Calculate significance using Fisher's exact tests
list.y.var <- unique(module_count_stat$Module)

Sub = list()
result = list()

for(i in list.y.var)
{
  sub = module_count_stat[module_count_stat$Module == i, c(5,6,2,3)]
  names(sub) <- rep(c("N","non_N"),2)
  sub = rbind(sub[,1:2],sub[,3:4])
  sub = as.matrix(sub)
  fisher <- fisher.test(sub)
  result[[i]] <- fisher.test(sub)
  Sub[[i]] <- fisher$p.value
}

fisher <- do.call("rbind", lapply(Sub, as.data.frame))
fisher <- data.frame(row.names(fisher),fisher)
names(fisher) <- c("Module", "Fisher_pval")

enrich_Macia_cycling <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_Macia_cycling$group <- "Macia_cycling"


# Look at genes with a time of day effect in Macia, that have no time of day effect in RTX (1468 genes in total)
#---------------------------------------------------------------------------------------------------------------
Macia_Time_DEGs <- dds_LRT_Time_Temperature.heat.Macia.signif
Macia_Time_DEGs$ID <- row.names(Macia_Time_DEGs)
Macia_Time_DEGs <- filter(Macia_Time_DEGs, ID %in% Macia_heat) #2690, ok
RTX_Time_DEGs <- dds_LRT_Time_Temperature.heat.RTX.signif
RTX_Time_DEGs$ID <- row.names(RTX_Time_DEGs)
RTX_Time_DEGs <- filter(RTX_Time_DEGs, ID %in% RTX_heat)
Macia_Time_DEGs_only <- anti_join(Macia_Time_DEGs, RTX_Time_DEGs, by = "ID")

setwd("C:/Users/tbonnot/Documents/sorghum/LRT_new/")
venn.diagram(list(B = Macia_Time_DEGs[,7], A = RTX_Time_DEGs[,7]), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_time_temp_genotypes.tiff")
write.table(Macia_Time_DEGs_only, "Macia_time_temp_only.txt", row.names = F, sep = "\t")

DEGs_to_look <- Macia_Time_DEGs_only$ID

DEGs_to_look_module <- filter(mod_Macia_heat, ID %in% DEGs_to_look)
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
Macia_module_count$Proportion_all <- (Macia_module_count$Freq_all)/(sum(Macia_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(Macia_module_count, DEGs_to_look_module_count, by = "Module")

# Calculate fold enrichment
module_count_stat$Fold_enrichment <- (module_count_stat$Proportion_DEG)/(module_count_stat$Proportion_all)

# Calculate significance using Fisher's exact tests
list.y.var <- unique(module_count_stat$Module)

Sub = list()
result = list()

for(i in list.y.var)
{
  sub = module_count_stat[module_count_stat$Module == i, c(5,6,2,3)]
  names(sub) <- rep(c("N","non_N"),2)
  sub = rbind(sub[,1:2],sub[,3:4])
  sub = as.matrix(sub)
  fisher <- fisher.test(sub)
  result[[i]] <- fisher.test(sub)
  Sub[[i]] <- fisher$p.value
}

fisher <- do.call("rbind", lapply(Sub, as.data.frame))
fisher <- data.frame(row.names(fisher),fisher)
names(fisher) <- c("Module", "Fisher_pval")

enrich_Macia_time_temp_only <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_Macia_time_temp_only$group <- "Macia_time_temp_only"


# Merge the different enrichment results
#---------------------------------------
enrich_DEG_results <- rbind(enrich_geno, enrich_Macia_TF,
                            enrich_Macia_cycling, enrich_Macia_time_temp_only)

# calculate the -log10(pval)
enrich_DEG_results$'|log10(Pval)|' <- -(log10(enrich_DEG_results$Fisher_pval))

# Indicate with Pval < 0.05
enrich_DEG_results <- enrich_DEG_results %>% mutate(Signif = ifelse(enrich_DEG_results$Fisher_pval < 0.05, "Yes", "No"))

# Define new module names
#------------------------
modules_names_Macia <- paste("Macia-H",seq(1,20,1), sep = "")
Module <- c("lightpink4","thistle1","brown","grey60", "saddlebrown", "darkorange", "coral1",
            "green", "mediumpurple3", "sienna3", "darkolivegreen", "cyan", 
            "plum", "black", "ivory", "darkslateblue", "honeydew1", "lightcyan1",
            "darkmagenta","greenyellow")
names_modules_Macia <- as.data.frame(cbind(modules_names_Macia, Module))
enrich_DEG_results <- left_join(enrich_DEG_results, names_modules_Macia, by = "Module")
enrich_DEG_results$modules_names_Macia <- factor(enrich_DEG_results$modules_names_Macia, levels = modules_names_Macia)

# plot the enrichment
#--------------------
enrich_DEG_results$group <- factor(enrich_DEG_results$group, levels = c("Geno","Macia_time_temp_only","Macia_TF","Macia_cycling"))
ggplot(enrich_DEG_results, 
       aes(x = modules_names_Macia, y = Fold_enrichment, group = group, colour = Signif)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=0.5)+
  geom_point(data=enrich_DEG_results,
             aes(x=modules_names_Macia, 
                 y=Fold_enrichment, group = group, fill = `|log10(Pval)|`, colour = Signif), 
             size = 4, shape = 21, stroke = 1)+
  scale_color_manual(values = c("grey","black"))+
  facet_wrap(~group, ncol = 1, scales = "free_y")+
  scale_fill_gradient(low="#FFECE6",high="#A42600",limits=c(0, 20), na.value="#420068ff")+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,0,5,0,"pt"), angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black", size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title.align=0.5)+
  xlab("GO biological processes")+
  ylab("Fold enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")


