###################################################
# Co-expression network analysis : cold stress data
###################################################
rm(list = ls())

library(WGCNA)
library(dplyr)

setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/cold_analysis/")

# Filter lowly expressed genes
#----------------------------
# Import raw counts
raw_counts <- read.table("countDFeByg.txt", header = T)

# Select cold stress data
raw_counts_cold <- raw_counts[,c(1:17,50:81)]

# keep genes with counts > 10 in at least 50% of the samples
raw_counts_cold$filter <- rowSums(raw_counts_cold[,-1] > 10,) >= 24
counts_cold_filtered <- raw_counts_cold[raw_counts_cold$filter == TRUE,]
#13780 genes

# Select rlog expression data for the filtered genes
#---------------------------------------------------
# Import rlog values
rlog <- read.table("rlog_expression_values.txt", header = T, stringsAsFactors = T)

# Select cold stress data
rlog2 <- rlog[,c(1:49)]

# Select filtered genes
rlog2 <- filter(rlog2, AGI %in% counts_cold_filtered$AGI)

# Separate by genotype 
RTX_cold <- rlog2[,c(1,grep("RTX", colnames(rlog2)))]
SC224_cold <- rlog2[,c(1,grep("SC224", colnames(rlog2)))]


# Network analysis by genotype
##############################
options(stringsAsFactors = FALSE)

# RTX430
########

RTXcoldData = RTX_cold

# Take a quick look at what is in the data set:
dim(RTXcoldData)
names(RTXcoldData)

# Remove the first columns corresponding to gene descriptions
datExpr0 = as.data.frame(t(RTXcoldData[, -1]))
names(datExpr0) = RTXcoldData$AGI
rownames(datExpr0) = names(RTXcoldData)[-1]

# Remove genes with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster the samples to see if there is any outlier
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(9,6)

pdf(file = "sampleClustering_RTX_cold.pdf", width = 9, height = 6);
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# No outlier

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# save the relevant expression and trait data for use in the next steps of the tutorial
save(datExpr, file = "RTX430-cold-dataInput.RData")

# Load the data saved in the first part
lnames = load(file = "RTX430-cold-dataInput.RData")

#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# Plot the results:
pdf("./scale_independence_RTX_cold.pdf", width=10, height=6)
#sizeGrWindow(12, 5)
par(mfrow = c(1,2))
cex1 = 0.6

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

# we choose the power 18 (see the vignette of WGCNA to know how to best choose this threshold)
softPower = 18
adjacency = adjacency(datExpr, power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merge modules that are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)

pdf(file = "geneDendro-RTX-cold.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "RTX430-cold-networkConstruction-stepByStep.RData")

# Load the expression and trait data saved in the first part
lnames = load(file = "RTX430-cold-dataInput.RData")

#The variable lnames contains the names of loaded variables.
lnames

# Load network data saved in the second part.
lnames = load(file = "RTX430-cold-networkConstruction-stepByStep.RData")
lnames

# Export the MEs info
write.table(MEs, "MEs_RTX_cold.txt", row.names = T, sep = "\t")

# import the annotation info
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/")
annot <- read.csv2("Sbicolor_454_v3.1.1.annotation_info_shiny.csv", header = T)
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/cold_analysis/")

# Select modules
colors_modules_all <- unique(moduleColors)
modules = colors_modules_all

# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = annot$Locus[match(modProbes, annot$Locus)]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               #edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               #nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


cyt_nodes <- cyt$nodeData
write.table(cyt_nodes, "modules_nodes_RTX_cold.txt", row.names = F, sep = "\t")

cyt_edges <- cyt$edgeData

# Export all edges with a weight > 0.15 
cyt_edges_0.15 <- cyt_edges[cyt_edges$weight > 0.15,]

names(cyt_nodes) <- c("fromNode","altName","module")
cyt_edges_0.15 <- merge.data.frame(cyt_edges_0.15[,1:3],cyt_nodes[,c(1,3)], by = "fromNode")
names(cyt_edges_0.15)[4] <- "module_from"
names(cyt_nodes) <- c("toNode","altName","module")
cyt_edges_0.15 <- merge.data.frame(cyt_edges_0.15, cyt_nodes[,c(1,3)], by = "toNode")

# Export the table to use in Cytoscape
write.table(cyt_edges_0.15, "RTX_cold_cytoscape_edges_0.15.txt", row.names = F, sep = "\t")

# Build a table with the module info and the gene names 
nodeAttr <- moduleColors[inModule]

# Check that the number of genes per module corresponds
RTX_mod_sum <- as.data.frame(table(nodeAttr)) # ok

# generate a table with the gene names and the module info
# careful, need to check after that because we are not merging, just assembling
# we guess that the order is the same, but it needs to be checked after that
RTX_mod <- data.frame(modProbes, moduleColors)

# export the module information
write.table(RTX_mod, "modules_RTX_cold.txt", row.names = F, sep = "\t")


# SC224
#######
rm(list = ls())
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/cold_analysis/")
SC224coldData = SC224_cold

# Take a quick look at what is in the data set:
dim(SC224coldData)
names(SC224coldData)

# Remove the first columns corresponding to gene descriptions
datExpr0 = as.data.frame(t(SC224coldData[, -1]))
names(datExpr0) = SC224coldData$AGI
rownames(datExpr0) = names(SC224coldData)[-1]

# Remove genes with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# If the last statement returns TRUE, all genes have passed the cuts
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster the samples to see if there is any outlier
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering_SC224_cold.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# No outlier

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# save the relevant expression and trait data for use in the next steps of the tutorial
save(datExpr0, file = "SC224-cold-dataInput.RData")

# Load the data saved in the first part
lnames = load(file = "SC224-cold-dataInput.RData")

#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, networkType = "signed")

# Plot the results:
pdf("./scale_independence_SC224_cold.pdf", width=10, height=6)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.6

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h = 0.8, col = "blue")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()


# we choose the power 18 (see the vignette of WGCNA to know how to best choose this threshold)
# We now calculate the adjacency, using the soft thresholding power 18:
softPower = 18
adjacency = adjacency(datExpr, power = softPower, type = "signed")

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merge modules that are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)

pdf(file = "geneDendro-SC224_cold.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "SC224-cold-networkConstruction-stepByStep.RData")

# Load the expression and trait data saved in the first part
lnames = load(file = "SC224-cold-dataInput.RData")

#The variable lnames contains the names of loaded variables.
lnames

# Load network data saved in the second part.
lnames = load(file = "SC224-cold-networkConstruction-stepByStep.RData")
lnames

# Export the MEs info
write.table(MEs, "MEs_SC224_cold.txt", row.names = T, sep = "\t")

# import the annotation info
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/")
annot <- read.csv2("Sbicolor_454_v3.1.1.annotation_info_shiny.csv", header = T)
setwd("C:/Users/tbonnot/Documents/sorghum/WGCNA/Signed_analysis/cold_analysis/")

# Select modules
colors_modules_all <- unique(moduleColors)
modules = colors_modules_all

# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = annot$Locus[match(modProbes, annot$Locus)]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               #edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               #nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

cyt_nodes <- cyt$nodeData
write.table(cyt_nodes, "modules_SC224_nodes_cold.txt", row.names = F, sep = "\t")
cyt_edges <- cyt$edgeData

# Export all edges with a weight > 0.15
cyt_edges_0.15 <- cyt_edges[cyt_edges$weight > 0.15,]

names(cyt_nodes) <- c("fromNode","altName","module")
cyt_edges_0.15 <- merge.data.frame(cyt_edges_0.15[,1:3],cyt_nodes[,c(1,3)], by = "fromNode")
names(cyt_edges_0.15)[4] <- "module_from"
names(cyt_nodes) <- c("toNode","altName","module")
cyt_edges_0.15 <- merge.data.frame(cyt_edges_0.15, cyt_nodes[,c(1,3)], by = "toNode")

# Export the table to use in Cytoscape
write.table(cyt_edges_0.15, "SC224_cold_cytoscape_edges_0.15.txt", row.names = F, sep = "\t")

# Build a table with the module info and the gene names (need to check after that!!)
nodeAttr <- moduleColors[inModule]

# Check that the number of genes per module corresponds
SC224_mod_sum <- as.data.frame(table(nodeAttr)) # looks good

# generate a table with the gene names and the module info
# careful, need to check after that because we are not merging, just assembling
# we guess that the order is the same, but it needs to be checked after that
SC224_mod <- data.frame(modProbes, moduleColors)

# export the module information
write.table(SC224_mod, "modules_SC224_cold.txt", row.names = F, sep = "\t")




