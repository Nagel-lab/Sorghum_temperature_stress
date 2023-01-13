############################################
# Analysis of DEGs from pairwise comparisons
############################################

library(ggplot2)
library(tidyr)
library(stringr)
library(VennDiagram) 
library(ade4)
library(factoextra)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(dplyr)
library(reshape2)

setwd("C:/Users/tbonnot/Documents/sorghum/Data_analysis/")

# import expression data and results of the pairwise comparisons
rlog <- read.table("rlog_expression_values.txt", header = T, stringsAsFactors = T)
rlog1 <- rlog[,c(1,50:length(rlog))]
rlog2 <- rlog[,c(1:49)]
LFC1 <- read.table("DESeq2_comp1.txt", header = T, stringsAsFactors = T)
LFC2 <- read.table("DESeq2_comp2.txt", header = T, stringsAsFactors = T)
LFC1$AGI <- row.names(LFC1)
LFC1 <- LFC1[,c(49,1:48)]
LFC2$AGI <- row.names(LFC2)
LFC2 <- LFC2[,c(49,1:48)]

# keep LFC and FDR columns in LFC datasets
data.LFC1 <- LFC1[,c(1,seq(3,45,by=6),seq(7,49,by=6))]
data.LFC1[is.na(data.LFC1)] <- 0
data.LFC2 <- LFC2[,c(1,seq(3,45,by=6),seq(7,49,by=6))]
data.LFC2[is.na(data.LFC2)] <- 0


# selection of DEGs (FDR < 5% and LFC >|1|)
#------------------------------------------
# Heat exp
H_RT.T1_high <- data.LFC1[(data.LFC1$H.RT.T1.H_H.RT.T1.C_FDR<0.05) & (data.LFC1$H.RT.T1.H_H.RT.T1.C_logFC>1),]
H_RT.T6_high <- data.LFC1[(data.LFC1$H.RT.T6.H_H.RT.T6.C_FDR<0.05) & (data.LFC1$H.RT.T6.H_H.RT.T6.C_logFC>1),]
H_RT.T9_high <- data.LFC1[(data.LFC1$H.RT.T9.H_H.RT.T9.C_FDR<0.05) & (data.LFC1$H.RT.T9.H_H.RT.T9.C_logFC>1),]
H_RT.T15_high <- data.LFC1[(data.LFC1$H.RT.T15.H_H.RT.T15.C_FDR<0.05) & (data.LFC1$H.RT.T15.H_H.RT.T15.C_logFC>1),]

H_RT.T1_low <- data.LFC1[(data.LFC1$H.RT.T1.H_H.RT.T1.C_FDR<0.05) & (data.LFC1$H.RT.T1.H_H.RT.T1.C_logFC<(-1)),]
H_RT.T6_low <- data.LFC1[(data.LFC1$H.RT.T6.H_H.RT.T6.C_FDR<0.05) & (data.LFC1$H.RT.T6.H_H.RT.T6.C_logFC<(-1)),]
H_RT.T9_low <- data.LFC1[(data.LFC1$H.RT.T9.H_H.RT.T9.C_FDR<0.05) & (data.LFC1$H.RT.T9.H_H.RT.T9.C_logFC<(-1)),]
H_RT.T15_low <- data.LFC1[(data.LFC1$H.RT.T15.H_H.RT.T15.C_FDR<0.05) & (data.LFC1$H.RT.T15.H_H.RT.T15.C_logFC<(-1)),]

H_MA.T1_high <- data.LFC1[(data.LFC1$H.MA.T1.H_H.MA.T1.C_FDR<0.05) & (data.LFC1$H.MA.T1.H_H.MA.T1.C_logFC>1),]
H_MA.T6_high <- data.LFC1[(data.LFC1$H.MA.T6.H_H.MA.T6.C_FDR<0.05) & (data.LFC1$H.MA.T6.H_H.MA.T6.C_logFC>1),]
H_MA.T9_high <- data.LFC1[(data.LFC1$H.MA.T9.H_H.MA.T9.C_FDR<0.05) & (data.LFC1$H.MA.T9.H_H.MA.T9.C_logFC>1),]
H_MA.T15_high <- data.LFC1[(data.LFC1$H.MA.T15.H_H.MA.T15.C_FDR<0.05) & (data.LFC1$H.MA.T15.H_H.MA.T15.C_logFC>1),]

H_MA.T1_low <- data.LFC1[(data.LFC1$H.MA.T1.H_H.MA.T1.C_FDR<0.05) & (data.LFC1$H.MA.T1.H_H.MA.T1.C_logFC<(-1)),]
H_MA.T6_low <- data.LFC1[(data.LFC1$H.MA.T6.H_H.MA.T6.C_FDR<0.05) & (data.LFC1$H.MA.T6.H_H.MA.T6.C_logFC<(-1)),]
H_MA.T9_low <- data.LFC1[(data.LFC1$H.MA.T9.H_H.MA.T9.C_FDR<0.05) & (data.LFC1$H.MA.T9.H_H.MA.T9.C_logFC<(-1)),]
H_MA.T15_low <- data.LFC1[(data.LFC1$H.MA.T15.H_H.MA.T15.C_FDR<0.05) & (data.LFC1$H.MA.T15.H_H.MA.T15.C_logFC<(-1)),]


# Cold exp
C_RT.T1_high <- data.LFC2[(data.LFC2$C.RT.T1.Cold_C.RT.T1.Cont_FDR<0.05) & (data.LFC2$C.RT.T1.Cold_C.RT.T1.Cont_logFC>1),]
C_RT.T6_high <- data.LFC2[(data.LFC2$C.RT.T6.Cold_C.RT.T6.Cont_FDR<0.05) & (data.LFC2$C.RT.T6.Cold_C.RT.T6.Cont_logFC>1),]
C_RT.T9_high <- data.LFC2[(data.LFC2$C.RT.T9.Cold_C.RT.T9.Cont_FDR<0.05) & (data.LFC2$C.RT.T9.Cold_C.RT.T9.Cont_logFC>1),]
C_RT.T15_high <- data.LFC2[(data.LFC2$C.RT.T15.Cold_C.RT.T15.Cont_FDR<0.05) & (data.LFC2$C.RT.T15.Cold_C.RT.T15.Cont_logFC>1),]

C_RT.T1_low <- data.LFC2[(data.LFC2$C.RT.T1.Cold_C.RT.T1.Cont_FDR<0.05) & (data.LFC2$C.RT.T1.Cold_C.RT.T1.Cont_logFC<(-1)),]
C_RT.T6_low <- data.LFC2[(data.LFC2$C.RT.T6.Cold_C.RT.T6.Cont_FDR<0.05) & (data.LFC2$C.RT.T6.Cold_C.RT.T6.Cont_logFC<(-1)),]
C_RT.T9_low <- data.LFC2[(data.LFC2$C.RT.T9.Cold_C.RT.T9.Cont_FDR<0.05) & (data.LFC2$C.RT.T9.Cold_C.RT.T9.Cont_logFC<(-1)),]
C_RT.T15_low <- data.LFC2[(data.LFC2$C.RT.T15.Cold_C.RT.T15.Cont_FDR<0.05) & (data.LFC2$C.RT.T15.Cold_C.RT.T15.Cont_logFC<(-1)),]

C_SC.T1_high <- data.LFC2[(data.LFC2$C.SC.T1.Cold_C.SC.T1.Cont_FDR<0.05) & (data.LFC2$C.SC.T1.Cold_C.SC.T1.Cont_logFC>1),]
C_SC.T6_high <- data.LFC2[(data.LFC2$C.SC.T6.Cold_C.SC.T6.Cont_FDR<0.05) & (data.LFC2$C.SC.T6.Cold_C.SC.T6.Cont_logFC>1),]
C_SC.T9_high <- data.LFC2[(data.LFC2$C.SC.T9.Cold_C.SC.T9.Cont_FDR<0.05) & (data.LFC2$C.SC.T9.Cold_C.SC.T9.Cont_logFC>1),]
C_SC.T15_high <- data.LFC2[(data.LFC2$C.SC.T15.Cold_C.SC.T15.Cont_FDR<0.05) & (data.LFC2$C.SC.T15.Cold_C.SC.T15.Cont_logFC>1),]

C_SC.T1_low <- data.LFC2[(data.LFC2$C.SC.T1.Cold_C.SC.T1.Cont_FDR<0.05) & (data.LFC2$C.SC.T1.Cold_C.SC.T1.Cont_logFC<(-1)),]
C_SC.T6_low <- data.LFC2[(data.LFC2$C.SC.T6.Cold_C.SC.T6.Cont_FDR<0.05) & (data.LFC2$C.SC.T6.Cold_C.SC.T6.Cont_logFC<(-1)),]
C_SC.T9_low <- data.LFC2[(data.LFC2$C.SC.T9.Cold_C.SC.T9.Cont_FDR<0.05) & (data.LFC2$C.SC.T9.Cold_C.SC.T9.Cont_logFC<(-1)),]
C_SC.T15_low <- data.LFC2[(data.LFC2$C.SC.T15.Cold_C.SC.T15.Cont_FDR<0.05) & (data.LFC2$C.SC.T15.Cold_C.SC.T15.Cont_logFC<(-1)),]

# count the number of DEGs by genotype
up_RTX_heat <- unique(c(H_RT.T1_high$AGI, H_RT.T6_high$AGI, H_RT.T9_high$AGI, H_RT.T15_high$AGI))
up_Macia_heat <- unique(c(H_MA.T1_high$AGI, H_MA.T6_high$AGI, H_MA.T9_high$AGI, H_MA.T15_high$AGI))
down_RTX_heat  <- unique(c(H_RT.T1_low$AGI, H_RT.T6_low$AGI, H_RT.T9_low$AGI, H_RT.T15_low$AGI))
down_Macia_heat <- unique(c(H_MA.T1_low$AGI, H_MA.T6_low$AGI, H_MA.T9_low$AGI, H_MA.T15_low$AGI))
up_RTX_cold <- unique(c(C_RT.T1_high$AGI, C_RT.T6_high$AGI, C_RT.T9_high$AGI, C_RT.T15_high$AGI))
up_SC224_cold <- unique(c(C_SC.T1_high$AGI, C_SC.T6_high$AGI, C_SC.T9_high$AGI, C_SC.T15_high$AGI))
down_RTX_cold <- unique(c(C_RT.T1_low$AGI, C_RT.T6_low$AGI, C_RT.T9_low$AGI, C_RT.T15_low$AGI))
down_SC224_cold <- unique(c(C_SC.T1_low$AGI, C_SC.T6_low$AGI, C_SC.T9_low$AGI, C_SC.T15_low$AGI))

# count the number of DEGs by time and genotype
comp <- c("H_RT.T1_high", "H_RT.T6_high", "H_RT.T9_high", "H_RT.T15_high",
          "H_RT.T1_low", "H_RT.T6_low", "H_RT.T9_low", "H_RT.T15_low",
          "H_MA.T1_high", "H_MA.T6_high", "H_MA.T9_high", "H_MA.T15_high",
          "H_MA.T1_low", "H_MA.T6_low", "H_MA.T9_low", "H_MA.T15_low",
          "C_RT.T1_high", "C_RT.T6_high", "C_RT.T9_high", "C_RT.T15_high",
          "C_RT.T1_low", "C_RT.T6_low", "C_RT.T9_low", "C_RT.T15_low",
          "C_SC.T1_high", "C_SC.T6_high", "C_SC.T9_high", "C_SC.T15_high",
          "C_SC.T1_low", "C_SC.T6_low", "C_SC.T9_low", "C_SC.T15_low")
count <- c(nrow(H_RT.T1_high), nrow(H_RT.T6_high), nrow(H_RT.T9_high), nrow(H_RT.T15_high),
           nrow(H_RT.T1_low), nrow(H_RT.T6_low), nrow(H_RT.T9_low), nrow(H_RT.T15_low),
           nrow(H_MA.T1_high), nrow(H_MA.T6_high), nrow(H_MA.T9_high), nrow(H_MA.T15_high),
           nrow(H_MA.T1_low), nrow(H_MA.T6_low), nrow(H_MA.T9_low), nrow(H_MA.T15_low),
           nrow(C_RT.T1_high), nrow(C_RT.T6_high), nrow(C_RT.T9_high), nrow(C_RT.T15_high),
           nrow(C_RT.T1_low), nrow(C_RT.T6_low), nrow(C_RT.T9_low), nrow(C_RT.T15_low),
           nrow(C_SC.T1_high), nrow(C_SC.T6_high), nrow(C_SC.T9_high), nrow(C_SC.T15_high),
           nrow(C_SC.T1_low), nrow(C_SC.T6_low), nrow(C_SC.T9_low), nrow(C_SC.T15_low))
count_DEGs <- data.frame(comp, count)

# Plot the number of DEGs by time of day and by genotype
#-------------------------------------------------------
count_DEGs <- count_DEGs %>% separate(comp, into = c("Exp","Genotype", "Time", "Response"))
count_DEGs$Time <- str_replace_all(count_DEGs$Time, 'T1', '1')
count_DEGs$Time <- str_replace_all(count_DEGs$Time, 'T6', '6')
count_DEGs$Time <- str_replace_all(count_DEGs$Time, 'T9', '9')
count_DEGs$Time <- str_replace_all(count_DEGs$Time, 'T15', '15')
count_DEGs$Genotype <- str_replace_all(count_DEGs$Genotype, 'MA', 'Macia')
count_DEGs$Genotype <- str_replace_all(count_DEGs$Genotype, 'RT', 'RTX430')
count_DEGs$Genotype <- str_replace_all(count_DEGs$Genotype, 'SC', 'SC224')
count_DEGs$Exp <- str_replace_all(count_DEGs$Exp, 'H', 'Heat experiment')
count_DEGs$Exp <- str_replace_all(count_DEGs$Exp, 'C', 'Cold experiment')
count_DEGs$Response <- str_replace_all(count_DEGs$Response, 'high', 'Up')
count_DEGs$Response <- str_replace_all(count_DEGs$Response, 'low', 'Down')

str(count_DEGs)
count_DEGs$Exp <- as.factor(count_DEGs$Exp)
count_DEGs$Genotype <- as.factor(count_DEGs$Genotype)
count_DEGs$Response <- as.factor(count_DEGs$Response)
count_DEGs$count <- as.numeric(as.character(count_DEGs$count))

# Heat experiment
count_DEGs_heat <- count_DEGs[count_DEGs$Exp == "Heat experiment",]
count_DEGs_heat$Time<- factor(count_DEGs_heat$Time,levels = c("1","6","9","15"))
count_DEGs_heat$Genotype<- factor(count_DEGs_heat$Genotype,levels = c("RTX430","Macia"))

ggplot(count_DEGs_heat, aes(x=Time)) +
  geom_bar(data=count_DEGs_heat[count_DEGs_heat$Response=="Up",], aes(y=count, fill= Response), stat="identity") +
  geom_bar(data=count_DEGs_heat[count_DEGs_heat$Response=="Down",], aes(y=-count, fill= Response), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  xlab("Time of Day")+
  ylab("Number of DEGs")+
  facet_wrap(~Genotype)+
  theme_bw()+
  scale_fill_manual(values = c("#10a9bcd1","#6e305cd1"))+#10bc6cd1
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        panel.grid = element_blank(),
        strip.background = element_blank())

# Cold experiment
count_DEGs_cold <- count_DEGs[count_DEGs$Exp == "Cold experiment",]
count_DEGs_cold$Time<- factor(count_DEGs_cold$Time,levels = c("1","6","9","15"))

ggplot(count_DEGs_cold, aes(x=Time)) +
  geom_bar(data=count_DEGs_cold[count_DEGs_cold$Response=="Up",], aes(y=count, fill= Response), stat="identity") +
  geom_bar(data=count_DEGs_cold[count_DEGs_cold$Response=="Down",], aes(y=-count, fill= Response), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  xlab("Time of Day")+
  ylab("Number of DEGs")+
  facet_wrap(~Genotype)+
  theme_bw()+
  scale_fill_manual(values = c("#10a9bcd1","#6e305cd1"))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        panel.grid = element_blank(),
        strip.background = element_blank())


# Export a table with the lists of all DEGs by time of day and genotype
#---------------------------------------------------------------------
L <- list(row.names(H_RT.T1_high), row.names(H_RT.T6_high), row.names(H_RT.T9_high), row.names(H_RT.T15_high),
          row.names(H_RT.T1_low), row.names(H_RT.T6_low), row.names(H_RT.T9_low), row.names(H_RT.T15_low),
          row.names(H_MA.T1_high), row.names(H_MA.T6_high), row.names(H_MA.T9_high), row.names(H_MA.T15_high),
          row.names(H_MA.T1_low), row.names(H_MA.T6_low), row.names(H_MA.T9_low), row.names(H_MA.T15_low),
          row.names(C_RT.T1_high), row.names(C_RT.T6_high), row.names(C_RT.T9_high), row.names(C_RT.T15_high),
          row.names(C_RT.T1_low), row.names(C_RT.T6_low), row.names(C_RT.T9_low), row.names(C_RT.T15_low),
          row.names(C_SC.T1_high), row.names(C_SC.T6_high), row.names(C_SC.T9_high), row.names(C_SC.T15_high),
          row.names(C_SC.T1_low), row.names(C_SC.T6_low), row.names(C_SC.T9_low), row.names(C_SC.T15_low))
cfun <- function(L) {
  pad.na <- function(x,len) {
    c(x,rep(NA,len-length(x)))
  }
  maxlen <- max(sapply(L,length))
  do.call(data.frame,lapply(L,pad.na,len=maxlen))
}
list_DEGs <- cfun(L)

names(list_DEGs) <- c("H_RT.T1_high", "H_RT.T6_high", "H_RT.T9_high", "H_RT.T15_high",
                      "H_RT.T1_low", "H_RT.T6_low", "H_RT.T9_low", "H_RT.T15_low",
                      "H_MA.T1_high", "H_MA.T6_high", "H_MA.T9_high", "H_MA.T15_high",
                      "H_MA.T1_low", "H_MA.T6_low", "H_MA.T9_low", "H_MA.T15_low",
                      "C_RT.T1_high", "C_RT.T6_high", "C_RT.T9_high", "C_RT.T15_high",
                      "C_RT.T1_low", "C_RT.T6_low", "C_RT.T9_low", "C_RT.T15_low",
                      "C_SC.T1_high", "C_SC.T6_high", "C_SC.T9_high", "C_SC.T15_high",
                      "C_SC.T1_low", "C_SC.T6_low", "C_SC.T9_low", "C_SC.T15_low")

setwd("C:/Users/tbonnot/Documents/sorghum/Data_analysis/")
write.table(list_DEGs, "list_DEGs.txt", row.names = F, sep = "\t")

# Export a table with the lists of DEGs by genotype
#--------------------------------------------------
RTX_heat_up <- list(row.names(H_RT.T1_high), row.names(H_RT.T6_high), row.names(H_RT.T9_high), row.names(H_RT.T15_high))
RTX_heat_up <- unique(unlist(RTX_heat_up))
RTX_heat_down <- list(row.names(H_RT.T1_low), row.names(H_RT.T6_low), row.names(H_RT.T9_low), row.names(H_RT.T15_low))
RTX_heat_down <- unique(unlist(RTX_heat_down))
Macia_heat_up <- list(row.names(H_MA.T1_high), row.names(H_MA.T6_high), row.names(H_MA.T9_high), row.names(H_MA.T15_high))
Macia_heat_up <- unique(unlist(Macia_heat_up))
Macia_heat_down <- list(row.names(H_MA.T1_low), row.names(H_MA.T6_low), row.names(H_MA.T9_low), row.names(H_MA.T15_low))
Macia_heat_down <- unique(unlist(Macia_heat_down))

RTX_cold_up <- list(row.names(C_RT.T1_high), row.names(C_RT.T6_high), row.names(C_RT.T9_high), row.names(C_RT.T15_high))
RTX_cold_up <- unique(unlist(RTX_cold_up))
RTX_cold_down <- list(row.names(C_RT.T1_low), row.names(C_RT.T6_low), row.names(C_RT.T9_low), row.names(C_RT.T15_low))
RTX_cold_down <- unique(unlist(RTX_cold_down))
SC224_cold_up <- list(row.names(C_SC.T1_high), row.names(C_SC.T6_high), row.names(C_SC.T9_high), row.names(C_SC.T15_high))
SC224_cold_up <- unique(unlist(SC224_cold_up))
SC224_cold_down <- list(row.names(C_SC.T1_low), row.names(C_SC.T6_low), row.names(C_SC.T9_low), row.names(C_SC.T15_low))
SC224_cold_down <- unique(unlist(SC224_cold_down))

L <- list(RTX_heat_up, RTX_heat_down, Macia_heat_up, Macia_heat_down, RTX_cold_up, RTX_cold_down, SC224_cold_up, SC224_cold_down)
cfun <- function(L) {
  pad.na <- function(x,len) {
    c(x,rep(NA,len-length(x)))
  }
  maxlen <- max(sapply(L,length))
  do.call(data.frame,lapply(L,pad.na,len=maxlen))
}
list_DEGs <- cfun(L)

names(list_DEGs) <- c("RTX_heat_up", "RTX_heat_down", "Macia_heat_up", "Macia_heat_down",
                      "RTX_cold_up", "RTX_cold_down", "SC224_cold_up", "SC224_cold_down")

setwd("C:/Users/tbonnot/Documents/sorghum/Data_analysis/")
write.table(list_DEGs, "DEGs_Fig1.txt", row.names = F, sep = "\t")


# Look at the overlapping DEGs between genotypes
#-----------------------------------------------
DEG_cold_RT_up <- rbind(C_RT.T1_high, C_RT.T6_high, C_RT.T9_high, C_RT.T15_high)
DEG_cold_RT_down <- rbind(C_RT.T1_low, C_RT.T6_low, C_RT.T9_low, C_RT.T15_low)
DEG_cold_SC_up <- rbind(C_SC.T1_high, C_SC.T6_high, C_SC.T9_high, C_SC.T15_high)
DEG_cold_SC_down <- rbind(C_SC.T1_low, C_SC.T6_low, C_SC.T9_low, C_SC.T15_low)

DEG_heat_RT_up <- rbind(H_RT.T1_high, H_RT.T6_high, H_RT.T9_high, H_RT.T15_high)
DEG_heat_RT_down <- rbind(H_RT.T1_low, H_RT.T6_low, H_RT.T9_low, H_RT.T15_low)
DEG_heat_Macia_up <- rbind(H_MA.T1_high, H_MA.T6_high, H_MA.T9_high, H_MA.T15_high)
DEG_heat_Macia_down <- rbind(H_MA.T1_low, H_MA.T6_low, H_MA.T9_low, H_MA.T15_low)

DEG_cold_RT_up <- DEG_cold_RT_up[!duplicated(DEG_cold_RT_up),]#735
DEG_cold_RT_down <- DEG_cold_RT_down[!duplicated(DEG_cold_RT_down),]#676
DEG_cold_SC_up <- DEG_cold_SC_up[!duplicated(DEG_cold_SC_up),]#963
DEG_cold_SC_down <- DEG_cold_SC_down[!duplicated(DEG_cold_SC_down),]#1125
DEG_heat_RT_up <- DEG_heat_RT_up[!duplicated(DEG_heat_RT_up),]#4024
DEG_heat_RT_down <- DEG_heat_RT_down[!duplicated(DEG_heat_RT_down),]#5177
DEG_heat_Macia_up <- DEG_heat_Macia_up[!duplicated(DEG_heat_Macia_up),]#3873
DEG_heat_Macia_down <- DEG_heat_Macia_down[!duplicated(DEG_heat_Macia_down),]#4197

overlap_cold_down <- merge.data.frame(DEG_cold_RT_down, DEG_cold_SC_down, by = "AGI")#349
overlap_cold_up <- merge.data.frame(DEG_cold_RT_up, DEG_cold_SC_up, by = "AGI")#491
overlap_heat_down <- merge.data.frame(DEG_heat_RT_down, DEG_heat_Macia_down, by = "AGI")#2816
overlap_heat_up <- merge.data.frame(DEG_heat_RT_up, DEG_heat_Macia_up, by = "AGI")#2900

venn.diagram(list(B = DEG_cold_RT_down[,1], A = DEG_cold_SC_down[,1]), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_cold_down.tiff")
venn.diagram(list(B = DEG_cold_RT_up[,1], A = DEG_cold_SC_up[,1]), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_cold_up.tiff")
venn.diagram(list(B = DEG_heat_RT_down[,1], A = DEG_heat_Macia_down[,1]), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_heat_down.tiff")
venn.diagram(list(B = DEG_heat_RT_up[,1], A = DEG_heat_Macia_up[,1]), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_heat_up.tiff")


# Filter expression data to perform PCA analyses from the identified DEGs
#########################################################################

DEG_heat <- rbind(H_RT.T1_high, H_RT.T6_high, H_RT.T9_high, H_RT.T15_high,
          H_RT.T1_low, H_RT.T6_low, H_RT.T9_low, H_RT.T15_low,
          H_MA.T1_high, H_MA.T6_high, H_MA.T9_high, H_MA.T15_high,
          H_MA.T1_low, H_MA.T6_low, H_MA.T9_low, H_MA.T15_low)
DEG_heat <- DEG_heat[!duplicated(DEG_heat),]
# 11,218 DEGs in the heat exp

DEG_cold <- rbind(C_RT.T1_high, C_RT.T6_high, C_RT.T9_high, C_RT.T15_high,
                  C_RT.T1_low, C_RT.T6_low, C_RT.T9_low, C_RT.T15_low,
                  C_SC.T1_high, C_SC.T6_high, C_SC.T9_high, C_SC.T15_high,
                  C_SC.T1_low, C_SC.T6_low, C_SC.T9_low, C_SC.T15_low)
DEG_cold <- DEG_cold[!duplicated(DEG_cold),]
# 2575 DEGs in the cold exp

# Select DEGs in the rlog table
#------------------------------
# Heat exp
rlog_heat_DEG <- merge.data.frame(rlog1, DEG_heat[,1:2], by = "AGI")
rlog_heat_DEG <- rlog_heat_DEG[,-50]

# Cold exp
rlog_cold_DEG <- merge.data.frame(rlog2, DEG_cold[,1:2], by = "AGI")
rlog_cold_DEG <- rlog_cold_DEG[,-50]

# Transpose the table to get the genes in columns
#------------------------------------------------
# Heat exp
rlog_heat_DEG_pca <- rlog_heat_DEG
row.names(rlog_heat_DEG_pca) <- rlog_heat_DEG_pca$AGI
rlog_heat_DEG_pca <- rlog_heat_DEG_pca[,-1]
rlog_heat_DEG_pca <- as.data.frame(t(rlog_heat_DEG_pca))
rlog_heat_DEG_pca$Condition <- row.names(rlog_heat_DEG_pca)
rlog_heat_DEG_pca <- rlog_heat_DEG_pca[,c(11219,1:11218)]

# Cold exp
rlog_cold_DEG_pca <- rlog_cold_DEG
row.names(rlog_cold_DEG_pca) <- rlog_cold_DEG_pca$AGI
rlog_cold_DEG_pca <- rlog_cold_DEG_pca[,-1]
rlog_cold_DEG_pca <- as.data.frame(t(rlog_cold_DEG_pca))
rlog_cold_DEG_pca$Condition <- row.names(rlog_cold_DEG_pca)
rlog_cold_DEG_pca <- rlog_cold_DEG_pca[,c(2576,1:2575)]

# Prepare the tables for PCA analyses
#------------------------------------
# heat
rlog_heat_DEG_pca <- rlog_heat_DEG_pca %>% separate(Condition, into = c("Exp","Genotype", "Time", "Temp"))
rlog_heat_DEG_pca <- rlog_heat_DEG_pca[,-1]
rlog_heat_DEG_pca$Genotype[rlog_heat_DEG_pca$Genotype == "RTX"] <- "RTX430"
rlog_heat_DEG_pca$Time <- str_replace_all(rlog_heat_DEG_pca$Time, 'T', '')
rlog_heat_DEG_pca$Temp <- str_replace_all(rlog_heat_DEG_pca$Temp, 'C1', '30.1')
rlog_heat_DEG_pca$Temp <- str_replace_all(rlog_heat_DEG_pca$Temp, 'C2', '30.2')
rlog_heat_DEG_pca$Temp <- str_replace_all(rlog_heat_DEG_pca$Temp, 'C3', '30.3')
rlog_heat_DEG_pca$Temp <- str_replace_all(rlog_heat_DEG_pca$Temp, 'H1', '42.1')
rlog_heat_DEG_pca$Temp <- str_replace_all(rlog_heat_DEG_pca$Temp, 'H2', '42.2')
rlog_heat_DEG_pca$Temp <- str_replace_all(rlog_heat_DEG_pca$Temp, 'H3', '42.3')
rlog_heat_DEG_pca <- rlog_heat_DEG_pca %>% separate(Temp, into = c("Temp","Rep"))

# cold
rlog_cold_DEG_pca <- rlog_cold_DEG_pca %>% separate(Condition, into = c("Exp","Genotype", "Time", "Temp"))
rlog_cold_DEG_pca <- rlog_cold_DEG_pca[,-1]
rlog_cold_DEG_pca$Genotype[rlog_cold_DEG_pca$Genotype == "RTX"] <- "RTX430"
rlog_cold_DEG_pca$Time <- str_replace_all(rlog_cold_DEG_pca$Time, 'T', '')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cont1', '30.1')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cont2', '30.2')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cont3', '30.3')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cold1', '10.1')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cold2', '10.2')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cold3', '10.3')
rlog_cold_DEG_pca$Temp <- str_replace_all(rlog_cold_DEG_pca$Temp, 'Cold4', '10.3')
rlog_cold_DEG_pca <- rlog_cold_DEG_pca %>% separate(Temp, into = c("Temp","Rep"))

# Perform the PCA
#----------------
## Heat exp
res.pca.heat <- dudi.pca(rlog_heat_DEG_pca[,-c(1:4)],nf=5, scannf = FALSE)
eig.val.heat <- get_eigenvalue(res.pca.heat)
head(eig.val.heat)

## Cold exp
res.pca.cold <- dudi.pca(rlog_cold_DEG_pca[,-c(1:4)],nf=5, scannf = FALSE)
eig.val.cold <- get_eigenvalue(res.pca.cold)
head(eig.val.cold)

# histogram to see the percentage of information per axis
#--------------------------------------------------------
fviz_screeplot(res.pca.heat,addlabels = TRUE, ncp=10, ylim=c(0,60),hjust = -0.1)
fviz_screeplot(res.pca.cold,addlabels = TRUE, ncp=10, ylim=c(0,60),hjust = -0.1)

# Visual representation of PCA
#------------------------------
## Heat exp
rlog_heat_DEG_pca$Time<- factor(rlog_heat_DEG_pca$Time,levels = c("1","6","9","15"))
rlog_heat_DEG_pca$Genotype<- factor(rlog_heat_DEG_pca$Genotype,levels = c("RTX430","Macia"))

# Axes 1 and 2
data_PCA_heat <- data.frame(res.pca.heat$li, rlog_heat_DEG_pca)
data_PCA_heat <- data_PCA_heat[,1:9]
data_PCA_heat$Time<- factor(data_PCA_heat$Time,levels = c("1","6","9","15"))

plot_pca_heat <- ggplot(data_PCA_heat, aes(x=Axis1, y=Axis2, color=Time, shape=interaction(Genotype,Temp), size = 0.5, stroke = 1.5))+
  geom_point()+
  theme_bw()+
  scale_shape_manual(values=c(21, 16, 24, 17))+
  scale_color_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18))

xdens <- axis_canvas(plot_pca_heat, axis = "x")+
  geom_density(data = data_PCA_heat, aes(x = Axis1, fill = Temp),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#b0b0b0", "#fd7700"))

ydens <- axis_canvas(plot_pca_heat, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_PCA_heat, aes(x = Axis2, fill = Time),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  scale_fill_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))

p1 <- insert_xaxis_grob(plot_pca_heat, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

# Axes 3 and 4
plot_pca_heat <- ggplot(data_PCA_heat, aes(x=Axis3, y=Axis4, color=Time, shape=interaction(Genotype,Temp), size = 0.5, stroke = 1.5))+
  geom_point()+
  theme_bw()+
  scale_shape_manual(values=c(21, 16, 24, 17))+
  scale_color_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18))

xdens <- axis_canvas(plot_pca_heat, axis = "x")+
  geom_density(data = data_PCA_heat, aes(x = Axis3, fill = Time),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))

ydens <- axis_canvas(plot_pca_heat, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_PCA_heat, aes(x = Axis4, fill = Genotype),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  scale_fill_manual(values = c("#d4d4d4",  "#3d3d3d"))

p1 <- insert_xaxis_grob(plot_pca_heat, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)


## Cold exp
rlog_cold_DEG_pca$Time<- factor(rlog_cold_DEG_pca$Time,levels = c("1","6","9","15"))

# Axes 1 and 2 
data_PCA_cold <- data.frame(res.pca.cold$li, rlog_cold_DEG_pca)
data_PCA_cold <- data_PCA_cold[,1:9]
data_PCA_cold$Time<- factor(data_PCA_cold$Time,levels = c("1","6","9","15"))

plot_pca_cold <- ggplot(data_PCA_cold, aes(x=Axis1, y=Axis2, color=Time, shape=interaction(Genotype,Temp), size = 0.5, stroke = 1.5))+
  geom_point()+
  theme_bw()+
  scale_shape_manual(values=c(23, 18, 21, 16))+
  scale_color_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18))

xdens <- axis_canvas(plot_pca_cold, axis = "x")+
  geom_density(data = data_PCA_cold, aes(x = Axis1, fill = Time),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))

ydens <- axis_canvas(plot_pca_cold, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_PCA_cold, aes(x = Axis2, fill = Temp),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  scale_fill_manual(values = c("#5693f9", "#b0b0b0"))

p1 <- insert_xaxis_grob(plot_pca_cold, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

# Axes 3 and 4 
plot_pca_cold <- ggplot(data_PCA_cold, aes(x=Axis3, y=Axis4, color=Time, shape=interaction(Genotype,Temp), size = 0.5, stroke = 1.5))+
  geom_point()+
  theme_bw()+
  scale_shape_manual(values=c(23, 18, 21, 16))+
  scale_color_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18))

xdens <- axis_canvas(plot_pca_cold, axis = "x")+
  geom_density(data = data_PCA_cold, aes(x = Axis3, fill = Time),
               alpha = 0.5, size = 0.2)+
  scale_fill_manual(values = c("#DDCC77","#88CCEE","#117733","#332288"))

ydens <- axis_canvas(plot_pca_cold, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_PCA_cold, aes(x = Axis4, fill = Genotype),
               alpha = 0.5, size = 0.2)+
  coord_flip()+
  scale_fill_manual(values = c("#d4d4d4", "#3d3d3d"))

p1 <- insert_xaxis_grob(plot_pca_cold, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)



# See how many genes are cycling within DEGs
############################################
Heat_up <- rbind(H_RT.T1_high, H_RT.T6_high, H_RT.T9_high, H_RT.T15_high,
                 H_MA.T1_high, H_MA.T6_high, H_MA.T9_high, H_MA.T15_high)
Heat_up <- Heat_up[!duplicated(Heat_up$AGI),1:2]
Heat_down <- rbind(H_RT.T1_low, H_RT.T6_low, H_RT.T9_low, H_RT.T15_low,
                 H_MA.T1_low, H_MA.T6_low, H_MA.T9_low, H_MA.T15_low)
Heat_down <- Heat_down[!duplicated(Heat_down$AGI),1:2]
Cold_up <- rbind(C_RT.T1_high, C_RT.T6_high, C_RT.T9_high, C_RT.T15_high,
                 C_SC.T1_high, C_SC.T6_high, C_SC.T9_high, C_SC.T15_high)
Cold_up <- Cold_up[!duplicated(Cold_up$AGI),1:2]
Cold_down <- rbind(C_RT.T1_low, C_RT.T6_low, C_RT.T9_low, C_RT.T15_low,
                   C_SC.T1_low, C_SC.T6_low, C_SC.T9_low, C_SC.T15_low)
Cold_down <- Cold_down[!duplicated(Cold_down$AGI),1:2]

# Import Lai et al data
setwd("C:/Users/tbonnot/Mes documents/sorghum/Phase_enrichment/")
Lai_cycling <- read.csv2("Phase_reference_Lai_rounded.csv", header = T, sep = ",")
names(Lai_cycling)[1] <- "AGI"
cycling_list <- Lai_cycling$AGI

Heat_up_cycling <- filter(Heat_up, AGI %in% cycling_list)
Heat_down_cycling <- filter(Heat_down, AGI %in% cycling_list)
Cold_up_cycling <- filter(Cold_up, AGI %in% cycling_list)
Cold_down_cycling <- filter(Cold_down, AGI %in% cycling_list)

# Prepare a table with the number of cycling gene by category of DEGs
Category <- c("Heat_up", "Heat_down", "Cold_up", "Cold_down")
DEGs <- c(nrow(Heat_up), nrow(Heat_down), nrow(Cold_up), nrow(Cold_down))
Cycling <- c(nrow(Heat_up_cycling), nrow(Heat_down_cycling), 
             nrow(Cold_up_cycling), nrow(Cold_down_cycling))
count_DEG_cycling <- as.data.frame(cbind(Category, DEGs, Cycling))
str(count_DEG_cycling)
count_DEG_cycling$DEGs <- as.numeric(count_DEG_cycling$DEGs)
count_DEG_cycling$Cycling <- as.numeric(count_DEG_cycling$Cycling)
count_DEG_cycling$Non_cycling <- (count_DEG_cycling$DEGs) - (count_DEG_cycling$Cycling)
library(reshape2)
count_DEG_cycling_plot <- melt(count_DEG_cycling[,-2], id.vars = "Category")

# Plot the proportion by group
count_DEG_cycling_plot$variable <- factor(count_DEG_cycling_plot$variable, levels = c("Non_cycling", "Cycling"))
ggplot(count_DEG_cycling_plot, aes(fill=variable, y=value, x=Category)) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        panel.border = element_blank())+
  scale_fill_manual(values = c("#c2c2c2","#117b1c"))


# See how many DEGs are TFs
#--------------------------
setwd("C:/Users/tbonnot/Documents/sorghum/PlantTFDB/")
TF_sorghum <- read.table("Sbi_TF_list.txt", header = T)

Cold_RTX430_up_TF <- filter(DEG_cold_RT_up, AGI %in% TF_sorghum$Gene_ID)
Cold_RTX430_down_TF <- filter(DEG_cold_RT_down, AGI %in% TF_sorghum$Gene_ID)
Cold_SC224_up_TF <- filter(DEG_cold_SC_up, AGI %in% TF_sorghum$Gene_ID)
Cold_SC224_down_TF <- filter(DEG_cold_SC_down, AGI %in% TF_sorghum$Gene_ID)

Heat_RTX430_up_TF <- filter(DEG_heat_RT_up, AGI %in% TF_sorghum$Gene_ID)
Heat_RTX430_down_TF <- filter(DEG_heat_RT_down, AGI %in% TF_sorghum$Gene_ID)
Heat_Macia_up_TF <- filter(DEG_heat_Macia_up, AGI %in% TF_sorghum$Gene_ID)
Heat_Macia_down_TF <- filter(DEG_heat_Macia_down, AGI %in% TF_sorghum$Gene_ID)

# Prepare a table with the number of TFs by category of DEGs
Category_2 <- c("Cold_RTX430_up", "Cold_RTX430_down","Cold_SC224_up","Cold_SC224_down",
                 "Heat_RTX430_up", "Heat_RTX430_down","Heat_Macia_up", "Heat_Macia_down")
DEGs_2 <- c(nrow(DEG_cold_RT_up), nrow(DEG_cold_RT_down), nrow(DEG_cold_SC_up), nrow(DEG_cold_SC_down),
            nrow(DEG_heat_RT_up), nrow(DEG_heat_RT_down), nrow(DEG_heat_Macia_up), nrow(DEG_heat_Macia_down))
DEG_TF <- c(nrow(Cold_RTX430_up_TF), nrow(Cold_RTX430_down_TF), 
            nrow(Cold_SC224_up_TF), nrow(Cold_SC224_down_TF),
            nrow(Heat_RTX430_up_TF), nrow(Heat_RTX430_down_TF), 
            nrow(Heat_Macia_up_TF), nrow(Heat_Macia_down_TF))
count_DEG_TF <- as.data.frame(cbind(Category_2, DEGs_2, DEG_TF))
str(count_DEG_TF)
count_DEG_TF$DEGs_2 <- as.numeric(count_DEG_TF$DEGs_2)
count_DEG_TF$DEG_TF <- as.numeric(count_DEG_TF$DEG_TF)
count_DEG_TF$Non_TF <- (count_DEG_TF$DEGs_2) - (count_DEG_TF$DEG_TF)

count_DEG_TF_plot <- melt(count_DEG_TF[,-2], id.vars = "Category_2")

# Plot the proportion by group
count_DEG_TF_plot$variable <- factor(count_DEG_TF_plot$variable, levels = c("Non_TF", "DEG_TF"))
ggplot(count_DEG_TF_plot, aes(fill=variable, y=value, x=Category_2)) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 12),
        panel.border = element_blank())+
  scale_fill_manual(values = c("#c2c2c2","#117b1c"))

