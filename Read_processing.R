##################################################################
# Read alignment, quantification and DEG analysis with SystemPipeR
##################################################################

# Please see the SystemPipeR vignette for details

library(systemPipeR)


## Read mapping with `HISAT2`
#----------------------------
args <- systemArgs(sysma="hisat2.param", mytargets="targets.txt")
sysargs(args)[1] # Command-line parameters for first FASTQ file
moduleload(modules(args))
system("hisat2-build ./Sbicolor_454_v3.fa ./Sbicolor_454_v3.fa")
resources <- list(walltime=120, ntasks=1, ncpus=5, memory=10240) 
reg <- clusterRun(args, conffile = ".batchtools.conf.R", Njobs=96, template = "batchtools.slurm.tmpl", runid="01", resourceList=resources)


## Read and alignment stats
#---------------------------
read_statsDF <- alignStats(args=args) 
write.table(read_statsDF, "alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


# Read quantification
#---------------------
## Read counting with `summarizeOverlaps` in parallel mode using multiple cores

library("GenomicFeatures"); library(BiocParallel)
txdb <- makeTxDbFromGFF(file="Sbicolor_454_v3.1.1.gene.gff3.gz", format="gff", dataSource="NCBI", organism="Sorghum bicolor")
saveDb(txdb, file="./sorbi.sqlite")
txdb <- loadDb("./sorbi.sqlite")
(align <- readGAlignments(outpaths(args)[1])) 
eByg <- exonsBy(txdb, by=c("gene"))
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=2); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=TRUE)) 
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
write.table(countDFeByg, "countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")


## Sample-wise correlation analysis
#----------------------------------
library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
countDF <- as.matrix(read.table("./countDFeByg.xls"))
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()


## Export rlog values
#--------------------
rlog <- assay(rlog(dds))
write.table(rlog, file = "rlog_expression_values.txt", sep = "\t")


## Analysis of DEGs (using DESeq2)
#-------------------
countDF <- read.delim("countDFeByg.xls", row.names=1, check.names=FALSE) 
targets <- read.delim("targets.txt", comment="#")
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
degseqDF1 <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE)
degseqDF2 <- run_DESeq2(countDF=countDF, targets=targets, cmp=cmp[[2]], independent=FALSE)
write.table(degseqDF1, "DESeq2_comp1.xls", quote=FALSE, sep="\t", col.names = NA)
write.table(degseqDF2, "DESeq2_comp2.xls", quote=FALSE, sep="\t", col.names = NA)


## Plot DEG results
#------------------
degseqDF1 <- read.delim("DESeq2_comp1.xls", row.names=1, check.names=FALSE) 
degseqDF2 <- read.delim("DESeq2_comp2.xls", row.names=1, check.names=FALSE) 

pdf("DEGcounts1.pdf")
DEG_list1 <- filterDEGs(degDF=degseqDF1, filter=c(Fold=2, FDR=5))
dev.off()
write.table(DEG_list1$Summary, "DEGcounts1.xls", quote=FALSE, sep="\t", row.names=FALSE)

pdf("DEGcounts2.pdf")
DEG_list2 <- filterDEGs(degDF=degseqDF2, filter=c(Fold=2, FDR=5))
dev.off()
write.table(DEG_list2$Summary, "DEGcounts2.xls", quote=FALSE, sep="\t", row.names=FALSE)


# Clustering and heat maps
#-------------------------
library(pheatmap)
geneids <- unique(as.character(unlist(DEG_list[[1]])))
y <- assay(rlog(dds))[geneids, ]
pdf("heatmap.pdf")
pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()

