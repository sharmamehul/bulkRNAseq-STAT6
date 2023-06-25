######################################
# 	 Load required libraries
######################################

library(pheatmap)
library(edgeR)
library(ggplot2)
library(limma)
library(data.table)
library(BiocManager)
library("org.Hs.eg.db")
library("AnnotationDbi")
library(pheatmap)


#convert ENSEMBL IDs to gene symbols
counts <- as.data.frame(counts)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts2 <- as.data.frame(counts)
counts2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(counts2), keytype = "ENSEMBL", column = "SYMBOL")

#remove NAs
counts3 <- counts2[!is.na(counts2$symbol),]

#aggregrate duplicates
counts3 %>% group_by(symbol) %>% summarise_all(sum) %>% data.frame() -> counts4
rownames(counts4) <- counts4$symbol
counts4 <- counts4[,-1]

#remove zeros
counts5<-data.matrix(counts4[,1:10])
dge  <-  DGEList(counts5, remove.zeros =   TRUE)
y<- dge[["counts"]]
y <- as.data.frame(y)
write.table(y, "cpmData.txt")

#load design matrix
des <- as.data.frame(design_)

#filter
filtergenes <- filterByExpr(y, des, min.count = 5, min.prop = 0.7)
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
filtered_norm_final <- filtered_norm[2:11]
write.table(filtered_norm_final, "filter_refined.txt")
write.table(filtered_norm_final, "filter_refined.csv")



mat <-  log2(filtered_norm_final + 1)

#batch <- c("A","B","C","A","B","C","A","B","C","A","B","C","A","B","C","A","B","C")

#batch_correct <- removeBatchEffect(counts, batch)

#write.table(batch_correct, "batch_corrected.txt")


#PCA calculation
matpca <- prcomp(t(mat), center=TRUE, scale. = TRUE)

pcaDF <- as.data.frame(matpca$x)
write.csv(pcaDF, "pcamatrix.csv")


#load data for limma
data <- read.table("C:/AMA/STAT6/RNAseq/HKseq/bulkseq/TPMs/filter_refined.csv")

meta <- as.data.frame(colnames(data))
meta$condition <- c("P_post1", "P_post", "P_post", "P_pos1", "P_post3", "HC", "HC", "HC", "HC", "HC")
meta$id <- c("pre_HC", "pre_HC1", "pre_HC1", "pre_HC1", "pre_HC2", "pre_HC", "pre_HC", "pre_HC", "pre_HC", "pre_HC")
rownames(meta) <- meta$`colnames(data)`
meta <- meta[, c("condition", "id")]

#patient vs HC comparison
test.met <- meta[meta$id == "pre_HC", ]
test.dat <- data[, rownames(test.met)]
design <- model.matrix(~ condition, data = test.met)
fit <- lmFit(test.dat, design)
fit <- eBayes(fit)
hits <- topTable(fit, coef = 2, p.value = 0.05, number = 10000)
write.csv(hits, "pre_HC_hits.csv")

all_genes <- rownames(hits)

#heatmap of limma hits
meta <- meta[1:10,]
#meta <- meta[c(16:18, 7:9, 10:12, 1:3, 13:15, 4:6), ]
#meta <- meta[c(15:17, 7:8), ]
data <- data[, rownames(meta)]
colorPalette <- c("blue", "blue", "white", "red", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(data[all_genes, ], annotation = meta[, -3], scale = "row", cluster_cols = F, fontsize_row = 5, color = colorPalette, border_color = "white")

