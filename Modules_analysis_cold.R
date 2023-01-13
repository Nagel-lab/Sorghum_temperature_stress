#######################
# Modules analysis cold
#######################

rm(list = ls())
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/cold_analysis/")

library(reshape2)
library(tidyr)
library(stringr)
library(doBy)
library(ggplot2)
library(dplyr)

# Plot the eigengenes for each module
#####################################
MEs_RTX_cold <- read.table("MEs_RTX_cold.txt", header = T)
MEs_SC224_cold <- read.table("MEs_SC224_cold.txt", header = T)
MEs_RTX_cold$sample <- row.names(MEs_RTX_cold)
MEs_SC224_cold$sample <- row.names(MEs_SC224_cold)

MEs_RTX_cold <- melt(MEs_RTX_cold, id.vars = "sample")
MEs_SC224_cold <- melt(MEs_SC224_cold, id.vars = "sample")
names(MEs_RTX_cold)[2] <- "module"
names(MEs_SC224_cold)[2] <- "module"
MEs_RTX_cold$Analysis <- "Modules_RTX"
MEs_SC224_cold$Analysis <- "Modules_SC224"

MEs_cold <- rbind(MEs_RTX_cold, MEs_SC224_cold)

# Separate the sample column
MEs_cold <- MEs_cold %>% separate(sample, into = c("Exp","Genotype", "Time", "Response"))
MEs_cold$Genotype <- str_replace_all(MEs_cold$Genotype, 'RTX', 'RTX430')
MEs_cold$Time <- str_replace_all(MEs_cold$Time, 'T', '')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cont1', '30.1')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cont2', '30.2')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cont3', '30.3')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cold1', '10.1')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cold2', '10.2')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cold3', '10.3')
MEs_cold$Response <- str_replace_all(MEs_cold$Response, 'Cold4', '10.4')
MEs_cold <- MEs_cold %>% separate(Response, into = c("Temperature", "Replicate"))

# calculate mean and SD by module
#--------------------------------
MEs_cold_mean <- summary_by(.~Analysis+Genotype+Time+Temperature+module, keep.names = T,FUN = function(x){mean(x,na.rm=T)} , data= MEs_cold)
MEs_cold_sd <- summary_by(.~Analysis+Genotype+Time+Temperature+module, keep.names = T,FUN = function(x){sd(x,na.rm=T)} , data= MEs_cold)
names(MEs_cold_sd)[6] <- "SD"
MEs_cold_mean <- merge.data.frame(MEs_cold_mean, MEs_cold_sd, by = c("Analysis","Genotype","Time","Temperature","module"))
names(MEs_cold_mean)[6] <- "Mean"

# Plot the eigengene by module and by analysis
#---------------------------------------------
str(MEs_cold_mean)
MEs_cold_mean$Analysis <- as.factor(MEs_cold_mean$Analysis)
MEs_cold_mean$Genotype <- factor(MEs_cold_mean$Genotype, levels = c("RTX430","SC224"))
MEs_cold_mean$Time <- as.numeric(as.character(MEs_cold_mean$Time))
MEs_cold_mean$Temperature <- as.factor(MEs_cold_mean$Temperature)

# Replace the default modules names by numbers from 1 to X (e.g. RTX-C1 for module 1 in RTX430 in the cold experiment)
Analysis_RTX <- MEs_cold_mean[MEs_cold_mean$Analysis == "Modules_RTX",]
modules_names_RTX <- paste("RTX-C",seq(1,19,1), sep = "")
module <- c("MEbrown","MEcoral2","MEcyan", "MEhoneydew1", "MEmediumpurple3", "MEdarkred",
                           "MElightsteelblue1", "MElightpink4", "MEdarkorange", "MEdarkorange2", "MEantiquewhite4", 
                           "MEgreenyellow", "MEthistle2", "MEred", "MEviolet", "MEcoral1", "MEgrey60", "MEdarkgreen", "MEgrey")
names_modules_RTX <- as.data.frame(cbind(modules_names_RTX, module))
Analysis_RTX <- left_join(Analysis_RTX, names_modules_RTX, by = "module")
Analysis_RTX$modules_names_RTX <- factor(Analysis_RTX$modules_names_RTX, levels = modules_names_RTX)

setwd("C:/Users/tbonnot/Documents/sorghum/Paper/Plots/")
pdf(file = "MEs_profiles_RTX_cold_figure.pdf", width = 10, height = 7)

# the module grey is removed because it corresponds to genes that were unsigned to any modules
ggplot(Analysis_RTX[!Analysis_RTX$module == "MEgrey",], aes(x = Time, y = Mean, group= Temperature)) +
  geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.08, color = NA)+
  geom_point(data=Analysis_RTX[!Analysis_RTX$module == "MEgrey",],aes(shape = Temperature, colour = Temperature),
             position=position_dodge(width=0.1), size = 3)+
  geom_line(data=Analysis_RTX[!Analysis_RTX$module == "MEgrey",],position=position_dodge(width=0.1),
            aes(linetype = Temperature, colour = Temperature))+
  scale_shape_manual(values=c(18,16))+
  scale_linetype_manual(values=c("blank", "solid"))+
  geom_errorbar(data=Analysis_RTX[!Analysis_RTX$module == "MEgrey",],
                aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                position=position_dodge(width=0.1))+
  scale_color_manual(values=c("#0066FF","#666666"))+
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

Analysis_SC224 <- MEs_cold_mean[MEs_cold_mean$Analysis == "Modules_SC224",]
modules_names_SC224 <- paste("SC224-C",seq(0,29,1), sep = "")
module <- c("MEgrey","MEcyan","MEskyblue2","MEskyblue3","MEplum2","MEsienna3","MEdarkturquoise",
            "MEpaleturquoise","MElightgreen", "MEfloralwhite","MEmaroon",
            "MEdarkgrey","MEdarkorange2","MEbisque4","MEwhite","MEdarkolivegreen",
            "MEmagenta","MElightcyan","MEthistle1","MEmediumpurple2",
            "MEgrey60","MEskyblue","MEmediumorchid","MEblack","MEorangered4","MEivory",
            "MEplum","MEdarkmagenta","MEcoral1","MEorangered3")

names_modules_SC224 <- as.data.frame(cbind(modules_names_SC224, module))
Analysis_SC224 <- left_join(Analysis_SC224, names_modules_SC224, by = "module")
Analysis_SC224$modules_names_SC224 <- factor(Analysis_SC224$modules_names_SC224, levels = modules_names_SC224)

pdf(file = "MEs_profiles_SC224_cold_figure.pdf", width = 10, height = 7)

ggplot(Analysis_SC224[!Analysis_SC224$module == "MEgrey",], aes(x = Time, y = Mean, group= Temperature)) +
  geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.08, color = NA)+
  geom_point(data=Analysis_SC224[!Analysis_SC224$module == "MEgrey",],aes(shape = Temperature, colour = Temperature),
             position=position_dodge(width=0.1), size = 3)+
  geom_line(data=Analysis_SC224[!Analysis_SC224$module == "MEgrey",],position=position_dodge(width=0.1),
            aes(linetype = Temperature, colour = Temperature))+
  scale_shape_manual(values=c(18,16))+
  scale_linetype_manual(values=c("blank", "solid"))+
  geom_errorbar(data=Analysis_SC224[!Analysis_SC224$module == "MEgrey",],aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                position=position_dodge(width=0.1))+
  scale_color_manual(values=c("#0066FF","#666666"))+
  scale_x_continuous(breaks = seq(0,15,3))+
  theme_bw()+
  facet_wrap(~modules_names_SC224)+
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
ggplot(Analysis_SC224[Analysis_SC224$module == "MEcyan",], aes(x = Time, y = Mean, group= Temperature)) +
  geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.08, color = NA)+
  geom_point(data=Analysis_SC224[Analysis_SC224$module == "MEcyan",],
             aes(shape = Temperature, colour = Temperature),
             position=position_dodge(width=0.1), size = 3)+
  geom_line(data=Analysis_SC224[Analysis_SC224$module == "MEcyan",],
            position=position_dodge(width=0.1),
            aes(linetype = Temperature, colour = Temperature))+
  scale_shape_manual(values=c(18, 16))+
  scale_linetype_manual(values=c("blank", "solid"))+
  geom_errorbar(data=Analysis_SC224[Analysis_SC224$module == "MEcyan",],
                aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                position=position_dodge(width=0.1))+
  scale_color_manual(values=c("#0066FF","#666666"))+
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
load(file = "dds_LRT_Temperature_Genotype_cold.RData")
load(file = "dds_LRT_Genotype.cold.RData")

setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/cold_analysis/")
dds_LRT_Temperature_Genotype.cold.signif <- dds_LRT_Temperature_Genotype.cold[dds_LRT_Temperature_Genotype.cold$padj < 0.05,]
dds_LRT_Temperature_Genotype.cold.signif <- na.omit(dds_LRT_Temperature_Genotype.cold.signif)
dds_LRT_Genotype.cold.signif <- dds_LRT_Genotype.cold[dds_LRT_Genotype.cold$padj < 0.05,]
dds_LRT_Genotype.cold.signif <- na.omit(dds_LRT_Genotype.cold.signif)

# count the number of genes by module in SC224
mod_RTX_cold <- read.table("modules_RTX_cold.txt", header = T)
mod_SC224_cold <- read.table("modules_SC224_cold.txt", header = T)

SC224_module_count <- data.frame(table(mod_SC224_cold$moduleColors))
names(SC224_module_count) <- c("Module","Freq_all")
SC224_module_count$Non_all <- sum(SC224_module_count$Freq_all)-(SC224_module_count$Freq_all)

# Look at the 5024 genes with a genotype effect and see if they are enriched in any module in SC224
#--------------------------------------------------------------------------------------------------
DEGs_to_look <- dds_LRT_Genotype.cold.signif
DEGs_to_look$ID <- row.names(DEGs_to_look)
names(mod_SC224_cold)[1] <- "ID"
DEGs_to_look_module <- merge.data.frame(mod_SC224_cold, DEGs_to_look, by = "ID")
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
SC224_module_count$Proportion_all <- (SC224_module_count$Freq_all)/(sum(SC224_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(SC224_module_count, DEGs_to_look_module_count, by = "Module")

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

DEGs_to_look_module <- filter(mod_SC224_cold, ID %in% DEGs_to_look)
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
SC224_module_count$Proportion_all <- (SC224_module_count$Freq_all)/(sum(SC224_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(SC224_module_count, DEGs_to_look_module_count, by = "Module")

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

enrich_SC224_TF <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_SC224_TF$group <- "SC224_TF"


# Look at cycling genes and see if they are enriched in any module in Macia
#--------------------------------------------------------------------------
setwd("C:/Users/tbonnot/Documents/sorghum/Phase_enrichment/")
Cycling_sorghum <- read.csv2("Phase_reference_Lai_rounded.csv", header = T, sep = ",")
DEGs_to_look <- Cycling_sorghum$ï..AGI

DEGs_to_look_module <- filter(mod_SC224_cold, ID %in% DEGs_to_look)
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion per module
SC224_module_count$Proportion_all <- (SC224_module_count$Freq_all)/(sum(SC224_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(SC224_module_count, DEGs_to_look_module_count, by = "Module")

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

enrich_SC224_cycling <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_SC224_cycling$group <- "SC224_cycling"


# Look at the genes with a time of day effect in SC224, that have no time of day effect in RTX (741 genes in total)
#------------------------------------------------------------------------------------------------------------------
# run LRT_analysis
SC224_Time_DEGs <- dds_LRT_Time_Temperature.cold.SC224.signif
SC224_Time_DEGs$ID <- row.names(SC224_Time_DEGs)
SC224_Time_DEGs <- filter(SC224_Time_DEGs, ID %in% SC224_cold)
RTX_Time_DEGs <- dds_LRT_Time_Temperature.cold.RTX.signif
RTX_Time_DEGs$ID <- row.names(RTX_Time_DEGs)
RTX_Time_DEGs <- filter(RTX_Time_DEGs, ID %in% RTX_cold)
SC224_Time_DEGs_only <- anti_join(SC224_Time_DEGs, RTX_Time_DEGs, by = "ID")

setwd("C:/Users/tbonnot/Documents/sorghum/LRT_new/")
venn.diagram(list(B = SC224_Time_DEGs[,7], A = RTX_Time_DEGs[,7]), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_time_temp_genotypes.tiff")
write.table(SC224_Time_DEGs_only, "SC224_time_temp_only.txt", row.names = F, sep = "\t")

DEGs_to_look <- SC224_Time_DEGs_only$ID


DEGs_to_look_module <- filter(mod_SC224_cold, ID %in% DEGs_to_look)
DEGs_to_look_module_count <- data.frame(table(DEGs_to_look_module$moduleColors))
names(DEGs_to_look_module_count) <- c("Module","Freq_DEG")
DEGs_to_look_module_count$Non_DEG <- sum(DEGs_to_look_module_count$Freq_DEG)-(DEGs_to_look_module_count$Freq_DEG)

# Calculate the proportion by module
SC224_module_count$Proportion_all <- (SC224_module_count$Freq_all)/(sum(SC224_module_count$Freq_all))
DEGs_to_look_module_count$Proportion_DEG <- (DEGs_to_look_module_count$Freq_DEG)/(sum(DEGs_to_look_module_count$Freq_DEG))

# Merge the reference with the DEG
module_count_stat <- merge.data.frame(SC224_module_count, DEGs_to_look_module_count, by = "Module")

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

enrich_SC224_time_temp_spe <- merge.data.frame(module_count_stat, fisher, by = "Module")
enrich_SC224_time_temp_spe$group <- "SC224_time_temp_spe"


# Merge the different enrichment results
#---------------------------------------
enrich_DEG_results <- rbind(enrich_geno, enrich_SC224_TF, enrich_SC224_cycling, enrich_SC224_time_temp_spe)

# calculate the -log10(pval)
enrich_DEG_results$'|log10(Pval)|' <- -(log10(enrich_DEG_results$Fisher_pval))

# Indicate with Pval < 0.05
enrich_DEG_results <- enrich_DEG_results %>% mutate(Signif = ifelse(enrich_DEG_results$Fisher_pval < 0.05, "Yes", "No"))

# Define new module names
#------------------------
modules_names_SC224 <- paste("SC224-C",seq(0,29,1), sep = "")
Module <- c("grey","cyan","skyblue2","skyblue3","plum2","sienna3","darkturquoise",
            "paleturquoise","lightgreen", "floralwhite","maroon",
            "darkgrey","darkorange2","bisque4","white","darkolivegreen",
            "magenta","lightcyan","thistle1","mediumpurple2",
            "grey60","skyblue","mediumorchid","black","orangered4","ivory",
            "plum","darkmagenta","coral1","orangered3")

names_modules_SC224_2 <- as.data.frame(cbind(modules_names_SC224, Module))
enrich_DEG_results <- left_join(enrich_DEG_results, names_modules_SC224_2, by = "Module")
enrich_DEG_results$modules_names_SC224 <- factor(enrich_DEG_results$modules_names_SC224, levels = modules_names_SC224)

# plot the enrichment
#--------------------
enrich_DEG_results$group <- factor(enrich_DEG_results$group, levels = c("Geno","SC224_time_temp_spe","SC224_TF","SC224_cycling"))
ggplot(enrich_DEG_results[!enrich_DEG_results$Module == "grey",], 
       aes(x = modules_names_SC224, y = Fold_enrichment, group = group, colour = Signif)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", size=0.5)+
  geom_point(data=enrich_DEG_results[!enrich_DEG_results$Module == "grey",],
             aes(x= modules_names_SC224, 
                 y=Fold_enrichment, group = group, fill = `|log10(Pval)|`, colour = Signif), 
             size = 4, shape = 21, stroke = 1)+
  scale_color_manual(values = c("grey","black"))+
  facet_wrap(~group, ncol = 1, scales = "free_y")+
  scale_fill_gradient(low="#EEF7FF",high="#2701A7",limits=c(0, 14), na.value="#420068ff")+
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

