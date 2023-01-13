###################
# LRT sorghum data
###################

setwd("C:/Users/tbonnot/Documents/sorghum/")

library(DESeq2)
library(tidyr)
library(stringr)

# Import raw counts
counts <- read.table("countDFeByg.txt", header = T)

# remove genes with 10 counts or less
counts$counts <- rowSums(counts[,2:length(counts)])
counts.clean <- counts[counts$counts >10,]
counts.clean <- counts.clean[,-97]
# 25532 genes with counts > 10

# Separate the count table into heat stress and cold stress experiments
counts.heat <- counts.clean[,c(17:48,81:96)]
counts.cold <- counts.clean[,c(1:16,49:80)]


# Prepare the metadata file
#--------------------------
# Heat
SampleName_heat <- names(counts.heat)
metadata.heat <- as.data.frame(SampleName_heat)
metadata.heat$SampleName <- metadata.heat$SampleName_heat
metadata.heat <- metadata.heat %>% separate(SampleName_heat, into = c("Exp","Genotype", "Time", "Response"))
metadata.heat$Response <- str_replace_all(metadata.heat$Response, 'C1', 'Cont.1')
metadata.heat$Response <- str_replace_all(metadata.heat$Response, 'C2', 'Cont.2')
metadata.heat$Response <- str_replace_all(metadata.heat$Response, 'C3', 'Cont.3')
metadata.heat$Response <- str_replace_all(metadata.heat$Response, 'H1', 'Heat.1')
metadata.heat$Response <- str_replace_all(metadata.heat$Response, 'H2', 'Heat.2')
metadata.heat$Response <- str_replace_all(metadata.heat$Response, 'H3', 'Heat.3')
metadata.heat <- metadata.heat %>% separate(Response, into = c("Temperature","Replicate"))
row.names(metadata.heat) <- metadata.heat$SampleName

# Cold
SampleName_cold <- names(counts.cold)
metadata.cold <- as.data.frame(SampleName_cold)
metadata.cold$SampleName <- metadata.cold$SampleName_cold
metadata.cold <- metadata.cold %>% separate(SampleName_cold, into = c("Exp","Genotype", "Time", "Response"))
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cont1', 'Cont.1')
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cont2', 'Cont.2')
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cont3', 'Cont.3')
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cold1', 'Cold.1')
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cold2', 'Cold.2')
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cold3', 'Cold.3')
metadata.cold$Response <- str_replace_all(metadata.cold$Response, 'Cold4', 'Cold.4')
metadata.cold <- metadata.cold %>% separate(Response, into = c("Temperature","Replicate"))
row.names(metadata.cold) <- metadata.cold$SampleName


# Check that the count tables and the metadata files are in the same order
#-------------------------------------------------------------------------
all(rownames(metadata.heat) %in% colnames(counts.heat))
all(rownames(metadata.heat) == colnames(counts.heat)) # TRUE

all(rownames(metadata.heat) %in% colnames(counts.heat))
all(rownames(metadata.heat) == colnames(counts.heat)) # TRUE


# LRT for the heat experiment
#############################
# Build the model
dds <- DESeqDataSetFromMatrix(countData = counts.heat,
                              colData = metadata.heat,
                              design = ~Time+Temperature+Genotype+Time:Temperature+Time:Genotype+Temperature:Genotype+Time:Temperature:Genotype)
# Test the effect of Time
dds_LRT_Time.heat <- DESeq(dds, test = "LRT", reduced = ~Temperature+Genotype+Temperature:Genotype)

# Test the effect of Genotype
dds_LRT_Genotype.heat <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature+Time:Temperature)

# Test the effect of Temperature
dds_LRT_Temperature.heat <- DESeq(dds, test = "LRT", reduced = ~Time+Genotype+Time:Genotype)

# Test the effect of Time:Temperature
dds_LRT_Time_Temperature.heat <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature+Genotype+Temperature:Genotype+Time:Genotype)

# Test the effect of Time:Genotype
dds_LRT_Time_Genotype.heat <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature+Genotype+Time:Temperature+Temperature:Genotype)

# Test the effect of Temperature:Genotype
dds_LRT_Temperature_Genotype.heat <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature+Genotype+Time:Temperature+Time:Genotype)

# Test the effect of Time:Temperature:Genotype
dds_LRT_Time_Temperature_Genotype.heat <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature+Genotype+Time:Temperature+Time:Genotype+Temperature:Genotype)


# Identify significant genes (FDR < 0.05) and exports the results
dds_LRT_Time.heat <- as.data.frame(results(dds_LRT_Time.heat))
dds_LRT_Time.heat.signif <- dds_LRT_Time.heat[dds_LRT_Time.heat$padj < 0.05,]
dds_LRT_Time.heat.signif <- na.omit(dds_LRT_Time.heat.signif)

dds_LRT_Genotype.heat <- as.data.frame(results(dds_LRT_Genotype.heat))
dds_LRT_Genotype.heat.signif <- dds_LRT_Genotype.heat[dds_LRT_Genotype.heat$padj < 0.05,]
dds_LRT_Genotype.heat.signif <- na.omit(dds_LRT_Genotype.heat.signif)

dds_LRT_Temperature.heat <- as.data.frame(results(dds_LRT_Temperature.heat))
dds_LRT_Temperature.heat.signif <- dds_LRT_Temperature.heat[dds_LRT_Temperature.heat$padj < 0.05,]
dds_LRT_Temperature.heat.signif <- na.omit(dds_LRT_Temperature.heat.signif)

dds_LRT_Time_Temperature.heat <- as.data.frame(results(dds_LRT_Time_Temperature.heat))
dds_LRT_Time_Temperature.heat.signif <- dds_LRT_Time_Temperature.heat[dds_LRT_Time_Temperature.heat$padj < 0.05,]
dds_LRT_Time_Temperature.heat.signif <- na.omit(dds_LRT_Time_Temperature.heat.signif)

dds_LRT_Time_Genotype.heat <- as.data.frame(results(dds_LRT_Time_Genotype.heat))
dds_LRT_Time_Genotype.heat.signif <- dds_LRT_Time_Genotype.heat[dds_LRT_Time_Genotype.heat$padj < 0.05,]
dds_LRT_Time_Genotype.heat.signif <- na.omit(dds_LRT_Time_Genotype.heat.signif)

dds_LRT_Temperature_Genotype.heat <- as.data.frame(results(dds_LRT_Temperature_Genotype.heat))
dds_LRT_Temperature_Genotype.heat.signif <- dds_LRT_Temperature_Genotype.heat[dds_LRT_Temperature_Genotype.heat$padj < 0.05,]
dds_LRT_Temperature_Genotype.heat.signif <- na.omit(dds_LRT_Temperature_Genotype.heat.signif)

dds_LRT_Time_Temperature_Genotype.heat <- as.data.frame(results(dds_LRT_Time_Temperature_Genotype.heat))
dds_LRT_Time_Temperature_Genotype.heat.signif <- dds_LRT_Time_Temperature_Genotype.heat[dds_LRT_Time_Temperature_Genotype.heat$padj < 0.05,]
dds_LRT_Time_Temperature_Genotype.heat.signif <- na.omit(dds_LRT_Time_Temperature_Genotype.heat.signif)

save(dds_LRT_Time.heat,
     file = "./LRT/dds_LRT_Time.heat.RData")
save(dds_LRT_Genotype.heat,
     file = "./LRT/dds_LRT_Genotype.heat.RData")
save(dds_LRT_Temperature.heat,
     file = "./LRT/dds_LRT_Temperature.heat.RData")
save(dds_LRT_Time_Temperature.heat,
     file = "./LRT/dds_LRT_Time_Temperature_heat.RData")
save(dds_LRT_Time_Genotype.heat,
     file = "./LRT/dds_LRT_Time_Genotype_heat.RData")
save(dds_LRT_Temperature_Genotype.heat,
     file = "./LRT/dds_LRT_Temperature_Genotype_heat.RData")
save(dds_LRT_Time_Temperature_Genotype.heat,
     file = "./LRT/dds_LRT_Time_Temperature_Genotype.heat.RData")


# LRT for the cold experiment
#############################
# Build the model
dds2 <- DESeqDataSetFromMatrix(countData = counts.cold,
                              colData = metadata.cold,
                              design = ~Time+Temperature+Genotype+Time:Temperature+Time:Genotype+Temperature:Genotype+Time:Temperature:Genotype)
# Test the effect of Time
dds_LRT_Time.cold <- DESeq(dds2, test = "LRT", reduced = ~Temperature+Genotype+Temperature:Genotype)

# Test the effect of Genotype
dds_LRT_Genotype.cold <- DESeq(dds2, test = "LRT", reduced = ~Time+Temperature+Time:Temperature)

# Test the effect of Temperature
dds_LRT_Temperature.cold <- DESeq(dds2, test = "LRT", reduced = ~Time+Genotype+Time:Genotype)

# Test the effect of Time:Temperature
dds_LRT_Time_Temperature.cold <- DESeq(dds2, test = "LRT", reduced = ~Time+Temperature+Genotype+Temperature:Genotype+Time:Genotype)

# Test the effect of Time:Genotype
dds_LRT_Time_Genotype.cold <- DESeq(dds2, test = "LRT", reduced = ~Time+Temperature+Genotype+Time:Temperature+Temperature:Genotype)

# Test the effect of Temperature:Genotype
dds_LRT_Temperature_Genotype.cold <- DESeq(dds2, test = "LRT", reduced = ~Time+Temperature+Genotype+Time:Temperature+Time:Genotype)

# Test the effect of Time:Temperature:Genotype
dds_LRT_Time_Temperature_Genotype.cold <- DESeq(dds2, test = "LRT", reduced = ~Time+Temperature+Genotype+Time:Temperature+Time:Genotype+Temperature:Genotype)


# Identify significant genes (FDR < 0.05) and exports the results
dds_LRT_Time.cold <- as.data.frame(results(dds_LRT_Time.cold))
dds_LRT_Time.cold.signif <- dds_LRT_Time.cold[dds_LRT_Time.cold$padj < 0.05,]
dds_LRT_Time.cold.signif <- na.omit(dds_LRT_Time.cold.signif)

dds_LRT_Genotype.cold <- as.data.frame(results(dds_LRT_Genotype.cold))
dds_LRT_Genotype.cold.signif <- dds_LRT_Genotype.cold[dds_LRT_Genotype.cold$padj < 0.05,]
dds_LRT_Genotype.cold.signif <- na.omit(dds_LRT_Genotype.cold.signif)

dds_LRT_Temperature.cold <- as.data.frame(results(dds_LRT_Temperature.cold))
dds_LRT_Temperature.cold.signif <- dds_LRT_Temperature.cold[dds_LRT_Temperature.cold$padj < 0.05,]
dds_LRT_Temperature.cold.signif <- na.omit(dds_LRT_Temperature.cold.signif)

dds_LRT_Time_Temperature.cold <- as.data.frame(results(dds_LRT_Time_Temperature.cold))
dds_LRT_Time_Temperature.cold.signif <- dds_LRT_Time_Temperature.cold[dds_LRT_Time_Temperature.cold$padj < 0.05,]
dds_LRT_Time_Temperature.cold.signif <- na.omit(dds_LRT_Time_Temperature.cold.signif)

dds_LRT_Time_Genotype.cold <- as.data.frame(results(dds_LRT_Time_Genotype.cold))
dds_LRT_Time_Genotype.cold.signif <- dds_LRT_Time_Genotype.cold[dds_LRT_Time_Genotype.cold$padj < 0.05,]
dds_LRT_Time_Genotype.cold.signif <- na.omit(dds_LRT_Time_Genotype.cold.signif)

dds_LRT_Temperature_Genotype.cold <- as.data.frame(results(dds_LRT_Temperature_Genotype.cold))
dds_LRT_Temperature_Genotype.cold.signif <- dds_LRT_Temperature_Genotype.cold[dds_LRT_Temperature_Genotype.cold$padj < 0.05,]
dds_LRT_Temperature_Genotype.cold.signif <- na.omit(dds_LRT_Temperature_Genotype.cold.signif)

dds_LRT_Time_Temperature_Genotype.cold <- as.data.frame(results(dds_LRT_Time_Temperature_Genotype.cold))
dds_LRT_Time_Temperature_Genotype.cold.signif <- dds_LRT_Time_Temperature_Genotype.cold[dds_LRT_Time_Temperature_Genotype.cold$padj < 0.05,]
dds_LRT_Time_Temperature_Genotype.cold.signif <- na.omit(dds_LRT_Time_Temperature_Genotype.cold.signif)

save(dds_LRT_Time.cold,
     file = "./LRT/dds_LRT_Time.cold.RData")
save(dds_LRT_Genotype.cold,
     file = "./LRT/dds_LRT_Genotype.cold.RData")
save(dds_LRT_Temperature.cold,
     file = "./LRT/dds_LRT_Temperature.cold.RData")
save(dds_LRT_Time_Temperature.cold,
     file = "./LRT/dds_LRT_Time_Temperature_cold.RData")
save(dds_LRT_Time_Genotype.cold,
     file = "./LRT/dds_LRT_Time_Genotype_cold.RData")
save(dds_LRT_Temperature_Genotype.cold,
     file = "./LRT/dds_LRT_Temperature_Genotype_cold.RData")
save(dds_LRT_Time_Temperature_Genotype.cold,
     file = "./LRT/dds_LRT_Time_Temperature_Genotype.cold.RData")


# LRT analyses to study the Time:Temperature interaction by genotype
####################################################################

# Heat experiment
#-----------------
metadata.heat.RTX <- metadata.heat[metadata.heat$Genotype == "RTX",]
counts.heat.RTX <- counts.heat[,grepl('RTX',names(counts.heat))]
metadata.heat.Macia <- metadata.heat[metadata.heat$Genotype == "Macia",]
counts.heat.Macia <- counts.heat[,grepl('Macia',names(counts.heat))]

all(rownames(metadata.heat.RTX) %in% colnames(counts.heat.RTX))
all(rownames(metadata.heat.RTX) == colnames(counts.heat.RTX)) # TRUE
all(rownames(metadata.heat.Macia) %in% colnames(counts.heat.Macia))
all(rownames(metadata.heat.Macia) == colnames(counts.heat.Macia)) # TRUE

# Build the model for RTX430 and test the effect of Time:Temperature
dds <- DESeqDataSetFromMatrix(countData = counts.heat.RTX,
                              colData = metadata.heat.RTX,
                              design = ~Time+Temperature+Time:Temperature)
dds_LRT_Time_Temperature.heat.RTX <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature)
dds_LRT_Time_Temperature.heat.RTX <- as.data.frame(results(dds_LRT_Time_Temperature.heat.RTX))

# Build the model for Macia and run the model
dds <- DESeqDataSetFromMatrix(countData = counts.heat.Macia,
                              colData = metadata.heat.Macia,
                              design = ~Time+Temperature+Time:Temperature)
dds_LRT_Time_Temperature.heat.Macia <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature)
dds_LRT_Time_Temperature.heat.Macia <- as.data.frame(results(dds_LRT_Time_Temperature.heat.Macia))


# Cold experiment
#-----------------
metadata.cold.RTX <- metadata.cold[metadata.cold$Genotype == "RTX",]
counts.cold.RTX <- counts.cold[,grepl('RTX',names(counts.cold))]
metadata.cold.SC224 <- metadata.cold[metadata.cold$Genotype == "SC224",]
counts.cold.SC224 <- counts.cold[,grepl('SC224',names(counts.cold))]

all(rownames(metadata.cold.RTX) %in% colnames(counts.cold.RTX))
all(rownames(metadata.cold.RTX) == colnames(counts.cold.RTX)) # TRUE
all(rownames(metadata.cold.SC224) %in% colnames(counts.cold.SC224))
all(rownames(metadata.cold.SC224) == colnames(counts.cold.SC224)) # TRUE

# Build the model for RTX430 and test the effect of Time:Temperature
dds <- DESeqDataSetFromMatrix(countData = counts.cold.RTX,
                              colData = metadata.cold.RTX,
                              design = ~Time+Temperature+Time:Temperature)

dds_LRT_Time_Temperature.cold.RTX <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature)
dds_LRT_Time_Temperature.cold.RTX <- as.data.frame(results(dds_LRT_Time_Temperature.cold.RTX))

# Build the model for SC224 and run the model
dds <- DESeqDataSetFromMatrix(countData = counts.cold.SC224,
                              colData = metadata.cold.SC224,
                              design = ~Time+Temperature+Time:Temperature)

dds_LRT_Time_Temperature.cold.SC224 <- DESeq(dds, test = "LRT", reduced = ~Time+Temperature)
dds_LRT_Time_Temperature.cold.SC224 <- as.data.frame(results(dds_LRT_Time_Temperature.cold.SC224))


# Identify significant genes (FDR < 0.05) and exports the results
#----------------------------------------------------------------
dds_LRT_Time_Temperature.heat.RTX.signif <- dds_LRT_Time_Temperature.heat.RTX[dds_LRT_Time_Temperature.heat.RTX$padj < 0.05,]
dds_LRT_Time_Temperature.heat.RTX.signif <- na.omit(dds_LRT_Time_Temperature.heat.RTX.signif)
dds_LRT_Time_Temperature.heat.Macia.signif <- dds_LRT_Time_Temperature.heat.Macia[dds_LRT_Time_Temperature.heat.Macia$padj < 0.05,]
dds_LRT_Time_Temperature.heat.Macia.signif <- na.omit(dds_LRT_Time_Temperature.heat.Macia.signif)
dds_LRT_Time_Temperature.cold.RTX.signif <- dds_LRT_Time_Temperature.cold.RTX[dds_LRT_Time_Temperature.cold.RTX$padj < 0.05,]
dds_LRT_Time_Temperature.cold.RTX.signif <- na.omit(dds_LRT_Time_Temperature.cold.RTX.signif)
dds_LRT_Time_Temperature.cold.SC224.signif <- dds_LRT_Time_Temperature.cold.SC224[dds_LRT_Time_Temperature.cold.SC224$padj < 0.05,]
dds_LRT_Time_Temperature.cold.SC224.signif <- na.omit(dds_LRT_Time_Temperature.cold.SC224.signif)

save(dds_LRT_Time_Temperature.heat.RTX,
     file = "./LRT/dds_LRT_Time.Temp.heat.RTX.RData")
save(dds_LRT_Time_Temperature.heat.Macia,
     file = "./LRT/dds_LRT_Time.Temp.heat.Macia.RData")
save(dds_LRT_Time_Temperature.cold.RTX,
     file = "./LRT/dds_LRT_Time.Temp.cold.RTX.RData")
save(dds_LRT_Time_Temperature.cold.SC224,
     file = "./LRT/dds_LRT_Time.Temp.cold.SC224.RData")


