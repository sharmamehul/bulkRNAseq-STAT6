######################################
# 	 Load required libraries
######################################

cat("loading ryesequire libraries...")
suppressPackageStartupMessages(library(package = "gdata"))
suppressPackageStartupMessages(library(package = "edgeR"))
suppressPackageStartupMessages(library(package = "ggplot2"))
suppressPackageStartupMessages(library(package = "pheatmap"))
suppressPackageStartupMessages(library(package = "grid"))
suppressPackageStartupMessages(library(package = "bioDist"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(package = "readxl"))
library(limma)
library(data.table)
library("scatterplot3d")
library(MUREN)

#merge files with raw counts
setwd("C:/AMA/STAT6/RNAseq/Jurkats/with_newWT2_stim_unstim")
nm <- list.files(path="C:/AMA/STAT6/RNAseq/Jurkats/with_newWT2_stim_unstim")
data <- do.call(cbind, lapply(nm, function(x) read.table(file=x)[, 2]))
colnames(data) <- basename(nm)
genes <- read.table(file ="Unstim-1144-1.counts.genes")[,1]
data <- cbind(genes, data)
counts_matrix <- as.data.frame(data)

rownames(counts_matrix) <- counts_matrix[,1]
counts_matrix <- counts_matrix[,-1]

n<-dim(counts_matrix)[1]
counts<-counts_matrix[1:(n-2),]

counts<-as.matrix(counts[,1:24])
write.csv(counts, "C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim/counts.csv")


setwd("C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim")
#re-load counts
counts <- as.data.frame(counts)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- as.matrix(counts[,1:24])

#Normalize using MUREN
counts_norm <- muren_norm(counts)
rownames(counts_norm) <- rownames(counts)
counts_norm <- as.data.frame(counts_norm)

# Normalized counts
dge  <-  DGEList(counts_norm, remove.zeros =   TRUE)
dge  <-  calcNormFactors(object = dge, method = "TMM")
normCounts  <- cpm(dge, normalized.lib.sizes=TRUE)

setwd("C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim")

# Write normlized counts matrix (count per millions)
y <- as.data.frame(normCounts)
write.table(normCounts, "cpmData.txt")

#load design matrix
des <- as.data.frame(design_)

filtergenes <- filterByExpr(normCounts, des, min.count = 5, min.prop = 0.7)
x <- as.data.frame(filtergenes)
write.table(filtergenes, "filter.txt")

genes <- rownames(y)
genes <- as.data.frame(genes)

a <- c(genes,x,y)
a <- as.data.frame(a)

filtered <- subset(a, filtergenes!="FALSE")
filtered_norm <- as.data.frame(filtered) 

rownames(filtered_norm) <- filtered_norm[,1]
filtered_norm <- filtered_norm[,-1]
filtered_norm_final <- filtered_norm[2:25]
write.table(filtered_norm_final, "filter_refined.txt")
write.table(filtered_norm_final, "filter_refined.csv")

data <- read.table("C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim/filter_refined.csv")

mat <-  log2(filtered_norm_final + 1)

#batch <- c("A","A","A","A","A","A","A","B","A","C","C","C","C","C","C","C","C","C")

#batch_correct <- removeBatchEffect(mat, batch)

#write.table(batch_correct, "batch_corrected.txt")

#mat <- as.data.frame(batch_correct)


#PCA calculation
matpca <- prcomp(t(mat), center=TRUE, scale. = TRUE)

pcaDF <- as.data.frame(matpca$x)

write.csv(pcaDF, "pcamatrix.csv")
write.csv(summary(matpca)$importance, 'overall_PCA_summary.csv')

pcaDF <- as.data.frame(pcamatrix)

colours <- c("dark blue", "purple", "dark green", "light blue", "violet", "green")
colours <- colours[pcaDF$Group]
shape <- c(17,15,16,17,15,16)
shape <- shape[pcaDF$Group]
scatterplot3d(pcaDF[,3:5], pch = shape, color = colours, cex.symbols = 1.5)
legend("bottom", legend = c("1144s", "1256s", "WTs", "1144-unstim", "1256-unstim", "WT-unstim"),
       col =  c("dark blue", "purple", "dark green", "light blue", "violet", "green"), pch = c(17,15,16,17,15,16), inset = -0.35, xpd = TRUE, horiz = TRUE)


## load data
data <- read.table("C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim/filter_refined.csv")

#data_wt <- as.data.frame(data[c(7:9, 16:18)])
#meta_wt <- as.data.frame(colnames(data_wt))
#meta_wt$id <- c("WT", "WT", "WT", "WT", "WT", "WT")
#meta_wt$stim <- c("stim", "stim", "stim", "unstim", "unstim", "unstim")
#meta_wt$day <- c("N1", "N3", "N1", "N2", "N3")
#rownames(meta_wt) <- meta_wt$`colnames(data_wt)`
#meta_wt <- meta_wt[, c("id", "stim")]

#test.met <- meta_wt[meta_wt$id == "WT", ]
#test.dat <- data_wt[, rownames(test.met)]
#design <- model.matrix(~ stim, data = test.met)
#fit <- lmFit(data_wt, design)
#fit <- eBayes(fit)
#WT_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)

meta <- as.data.frame(colnames(data))
meta$stim <- c("stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "unstim")
meta$id <- c("1144", "1144", "1144", "1256", "1256", "1256", "1784", "1784", "1784", "WT", "WT", "WT", "1144", "1144", "1144", "1256", "1256", "1256", "1784", "1784", "1784", "WT", "WT", "WT")
#meta$day <- c("N1", "N2", "N3", "N1", "N2", "N3", "N1", "N3", "N1", "N2", "N3", "N1", "N2", "N3", "N1", "N2", "N3")
rownames(meta) <- meta$`colnames(data)`
meta <- meta[, c("stim", "id")]

#stim vs unstim comparison
test.met <- meta[meta$id == "1256", ]
test.dat <- data[, rownames(test.met)]
design <- model.matrix(~ stim, data = test.met)
fit <- lmFit(test.dat, design)
fit <- eBayes(fit)
hits_1256 <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)
write.csv(hits_1256, "1256_stimvsunstim.csv")

test.met <- meta[meta$id == "1144", ]
test.dat <- data[, rownames(test.met)]
design <- model.matrix(~ stim, data = test.met)
fit <- lmFit(test.dat, design)
fit <- eBayes(fit)
hits_1144 <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)
write.csv(hits_1144, "1144_stimvsunstim.csv")

test.met <- meta[meta$id == "1784", ]
test.dat <- data[, rownames(test.met)]
design <- model.matrix(~ stim, data = test.met)
fit <- lmFit(test.dat, design)
fit <- eBayes(fit)
hits_1784 <- topTable(fit, coef = 2, p.value = 0.2, number = 10000)
write.csv(hits_1784, "1784_stimvsunstim.csv")

test.met <- meta[meta$id == "WT", ]
test.dat <- data[, rownames(test.met)]
design <- model.matrix(~ stim, data = test.met)
fit <- lmFit(test.dat, design)
fit <- eBayes(fit)
WT_hits <- topTable(fit, coef = 2, p.value = 0.1, number = 10000)

#x <- as.data.frame(WT_hits[1]>0)
#b <- c(WT_hits, x)
#b <- as.data.frame(b)
#rownames(b) <- rownames(WT_hits)
#WT_up <- subset(b, logFC.1!="FALSE")
#WT_up <- as.data.frame(WT_up) 

write.csv(WT_hits, "WT_stimvsunstim_q0_1.csv")

#all_genes <- c(rownames(WT_hits), rownames(hits_1144), rownames(hits_1256), rownames(hits_1784)) %>% unique()
all_genes <- rownames(WT_hits)
#write.csv(WT_hits, "WT_stimvsunstim.csv")


#test.met <- meta[meta$stim == "stim", ]
#test.dat <- data[, rownames(test.met)]
#design <- model.matrix(~ id + day, data = test.met)
#fit <- lmFit(test.dat, design)
#fit <- eBayes(fit)
#var_hits <- topTable(fit, coef = 2, p.value = 0.2, number = 10000)

#all_genes <- rownames(var_hits)

#import gene list to use for heatmap
genes <- as.data.frame(Genes_stimvsunstim_adjpval05_fc125)
#genes <- as.data.frame(test_genes)
all_genes <- genes[,1]



library(pheatmap)
meta <- meta[c(22:24, 10:12, 13:15, 1:3, 16:18, 4:6, 19:21, 7:9),]
#meta <- meta[c(16:18, 7:9, 10:12, 1:3, 13:15, 4:6), ]
#meta <- meta[c(15:17, 7:8), ]
data <- data[, rownames(meta)]
colorPalette <- c("blue", "blue", "white", "red", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(data[all_genes, ], annotation = meta[, -3], scale = "row", cluster_cols = F, fontsize_row = 7, color = colorPalette, border_color = "white")








## load data _ for unsti comp
data <- read.table("C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim_no1784/filter_refined.txt")
data1 <- as.data.frame(colnames(data))

wt_1784_unstim <- as.data.frame(data[c(7:9, 10:12)])
meta_wt <- as.data.frame(colnames(wt_1784_unstim))
meta_wt$id <- c("WT", "WT", "WT", "1784", "1784", "1784")
meta_wt$stim <- c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1784_unstim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "unstim", ]
test.dat <- wt_1784_unstim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1784_unstim, design)
fit <- eBayes(fit)
wt_1784_unstim_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)
write.csv(wt_1784_unstim_hits, "wt_1784_unstim_hits.csv")


#wt vs 1256 unstim
wt_1256_unstim <- as.data.frame(data[c(13:15, 16:18)])
meta_wt <- as.data.frame(colnames(wt_1256_unstim))
meta_wt$id <- c("1256", "1256", "1256", "WT", "WT", "WT")
meta_wt$stim <- c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1256_unstim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "unstim", ]
test.dat <- wt_1256_unstim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1256_unstim, design)
fit <- eBayes(fit)
wt_1256_unstim_hits <- topTable(fit, coef = 2, p.value = 0.5, number = 10000)

write.csv(wt_1256_unstim_hits, "wt_1256_unstim.csv")


#wt vs 1144 unstim
wt_1144_unstim <- as.data.frame(data[c(10:12, 16:18)])
meta_wt <- as.data.frame(colnames(wt_1144_unstim))
meta_wt$id <- c("1144", "1144", "1144", "WT", "WT", "WT")
meta_wt$stim <- c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1144_unstim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "unstim", ]
test.dat <- wt_1144_unstim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1144_unstim, design)
fit <- eBayes(fit)
wt_1144_unstim_hits <- topTable(fit, coef = 2, p.value = 0.5, number = 10000)

write.csv(wt_1144_unstim_hits, "wt_1144_unstim.csv")


#wt vs 1144 unstim
wt_1144_stim <- as.data.frame(data[c(1:3, 7:9)])
meta_wt <- as.data.frame(colnames(wt_1144_stim))
meta_wt$id <- c("1144", "1144", "1144", "WT", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1144_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_1144_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1144_stim, design)
fit <- eBayes(fit)
wt_1144_stim_hits <- topTable(fit, coef = 2, p.value = 0.5, number = 10000)

write.csv(wt_1144_stim_hits, "wt_1144_stim.csv")

#wt vs 1256 unstim
wt_1256_stim <- as.data.frame(data[c(4:6, 7:9)])
meta_wt <- as.data.frame(colnames(wt_1256_stim))
meta_wt$id <- c("1256", "1256", "1256", "WT", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1256_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_1256_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1256_stim, design)
fit <- eBayes(fit)
wt_1256_stim_hits <- topTable(fit, coef = 2, p.value = 0.5, number = 10000)

write.csv(wt_1256_stim_hits, "wt_1256_stim.csv")


all_genes <- c(rownames(wt_1144_unstim_hits), rownames(wt_1256_unstim_hits)) %>% unique()

library(pheatmap)
data1 <- data1[c(15:17, 9:11, 12:14),]
data1 <- as.data.frame(data1)
rownames(data1) <- data1[,1]
data <- data[, rownames(data1)]
colorPalette <- c("red", "red", "white", "blue", "blue")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(data[all_genes, ], annotation = meta[, -3], scale = "row", cluster_cols = F, fontsize_row = 12, color = colorPalette, border_color = "white")




## load data _ for stim comp
data <- read.table("C:/AMA/STAT6/RNAseq/Jurkats/analysis_stim_unstim_without WT2sand1784/filter_refined.txt")
data1 <- as.data.frame(colnames(data))

wt_1144_stim <- as.data.frame(data[c(1:3, 7:8)])
meta_wt <- as.data.frame(colnames(wt_1144_stim))
meta_wt$id <- c("1144", "1144", "1144", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1144_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_1144_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1144_stim, design)
fit <- eBayes(fit)
wt_1144_stim_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)
write.csv(wt_1144_stim_hits, "wt_1144_stim_hits.csv")

wt_1256_stim <- as.data.frame(data[c(4:6, 7:8)])
meta_wt <- as.data.frame(colnames(wt_1256_stim))
meta_wt$id <- c("1256", "1256", "1256", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1256_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_1256_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1256_stim, design)
fit <- eBayes(fit)
wt_1256_stim_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)
write.csv(wt_1256_stim_hits, "wt_1256_stim_hits.csv")

all_genes <- c(rownames(wt_1144_stim_hits), rownames(wt_1256_stim_hits)) %>% unique()

library(pheatmap)
data1 <- data1[c(15:17, 7:8, 9:11, 1:3, 12:14, 4:6),]
data1 <- as.data.frame(data1)
rownames(data1) <- data1[,1]
data <- data[, rownames(data1)]
colorPalette <- c("red", "red", "white", "blue", "blue")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(data[all_genes, ], annotation = meta[, -3], scale = "row", cluster_cols = F, fontsize_row = 12, color = colorPalette, border_color = "white")






## load data _ compare everything to WT unstim
data <- read.table("C:/AMA/STAT6/RNAseq/Jurkats/analysis_stim_unstim_without WT2sand1784/filter_refined.txt")
data <- read.table("C:/AMA/STAT6/RNAseq/Jurkats/analysis_stim_unstim_without WT2sand1784/avg_counts.txt")
data1 <- as.data.frame(colnames(data))

#1144s
wt_1144_stim <- as.data.frame(data[c(1:3, 9:11)])
meta_wt <- as.data.frame(colnames(wt_1144_stim))
meta_wt$id <- c("1144", "1144", "1144", "WT", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1144_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_1144_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1144_stim, design)
fit <- eBayes(fit)
wt_1144_stim_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 100)

#1256s
wt_1256_stim <- as.data.frame(data[c(4:6, 9:11)])
meta_wt <- as.data.frame(colnames(wt_1256_stim))
meta_wt$id <- c("1256", "1256", "1256", "WT", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_1256_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_1256_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_1256_stim, design)
fit <- eBayes(fit)
wt_1256_stim_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 100)

#WTs
wt_WT_stim <- as.data.frame(data[c(7:8, 9:11)])
meta_wt <- as.data.frame(colnames(wt_WT_stim))
meta_wt$id <- c("WTs", "WTs", "WT", "WT", "WT")
meta_wt$stim <- c("stim", "stim", "stim", "stim", "stim")
rownames(meta_wt) <- meta_wt$`colnames(wt_WT_stim)`
meta_wt <- meta_wt[, c("id", "stim")]

test.met <- meta_wt[meta_wt$stim == "stim", ]
test.dat <- wt_WT_stim[, rownames(test.met)]
design <- model.matrix(~ id, data = test.met)
fit <- lmFit(wt_WT_stim, design)
fit <- eBayes(fit)
wt_WT_stim_hits <- topTable(fit, coef = 2, p.value = 0.05, number = 100)

all_genes <- c(rownames(wt_1144_stim_hits), rownames(wt_1256_stim_hits), rownames(wt_WT_stim_hits)) %>% unique()
all_genes <- c(rownames(data))


library(pheatmap)
data1 <- data1[c(6, 3, 4, 1, 5, 2),]
data1 <- as.data.frame(data1)
rownames(data1) <- data1[,1]
data <- data[, rownames(data1)]
colorPalette <- c("red", "red", "white", "blue", "blue")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(data[all_genes, ], scale = "row", cluster_cols = F, clustering_method	= "ward", fontsize_row = 12, color = colorPalette, border_color = "white", )






## SLEA function
# do SLEA (Sample Level Enrichment Analysis) function(z-score per sample)
doSLEA <- function(expressionSet, geneSet) {
  # scale expression
  exprsMat <- expressionSet
  exprsMat <- t(scale(t(exprsMat)))
  # extract expression of leGenes of each geneset
  comm <- intersect(geneSet, rownames(expressionSet))
  gsDF <- exprsMat[comm, ]
  # calculate mean expression per sample
  gsM <- colMeans(gsDF)
  # extract random genes of size of the geneSet from full probeset and calculate mean
  # and perform this for 'n' permutations
  nperm <- lapply(1:1000, function(j) {
    # set seed for every permutation
    set.seed(j)
    rGSDF <- exprsMat[sample.int(nrow(exprsMat),length(comm)), ]
    rGSM <- colMeans(rGSDF)
    return(value = rGSM)
  })
  permDF <- do.call(rbind, nperm)
  zscore <- (gsM - colMeans(permDF)) / apply(permDF,2,sd)
  sleaDF <- zscore %>% as.data.frame()
  return(value = sleaDF)
}

FILE = "C:/AMA/STAT6/RNAseq/Jurkats/analysis_with_newWT2_stim_unstim_no1784/gsea/pathways_hallmarks.xlsx"
gs1 <- read_excel(FILE, sheet = 1) %>%
  as.data.frame()
le1 <- as.list(gs1[,1:6])

# load expression file
the_file <- as.data.frame(counts_forSLEA)
rownames(the_file) <- as.character(unlist(the_file[1]))
the_file = the_file[,-1]
mat <- the_file

#if row mean is zero, remove row, and remove gene name for that row from genesets


# call SLEA
sleaLS <- lapply(1:length(le1), function(l) {
  expressionSet = mat
  geneSet <- le1[[l]] %>% strsplit(",") %>% unlist(.)
  sDF <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
  names(sDF) <- names(le1[l])
  return(value = sDF)
})
sleaDF <- do.call(cbind, sleaLS)
sleaDF <- sleaDF %>% t() %>% as.data.frame()

colorPalette <- c("blue", "white", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)

pheatmap(mat = sleaDF,
         color = colorPalette,
         cellwidth = 25,
         cellheight = 25,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = 0,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 10,
         fontsize_row = 14,
         border_color = "white",
         gaps_col = c(9))
dev.off()