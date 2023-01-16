## Loading libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(reticulate)
library(ggridges)
library(reshape2)

# Loading and pre-processing KalwRij SS 1 nk/t
K_ss_1_nkt <- Read10X(data.dir = "~/kal_1_ss/filtered_feature_bc_matrix/")
K_ss_1_nkt <- CreateSeuratObject(counts = K_ss_1_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=K_ss_1_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = K_ss_1_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = K_ss_1_nkt, slot = "counts"))
K_ss_1_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = K_ss_1_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
K_ss_1_nkt <- subset(x = K_ss_1_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_K_ss_1_nkt <- as.character(read.csv(file = "~/barcodes_kss1_nkt.csv", header = T, row.names = 1)[,1])
nk_K_ss_1_nkt <- gsub("_.", "", nk_K_ss_1_nkt)
K_ss_1_nkt <- subset(x=K_ss_1_nkt, cells = nk_K_ss_1_nkt)
K_ss_1_nkt <- NormalizeData(object = K_ss_1_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
K_ss_1_nkt <- FindVariableFeatures(object = K_ss_1_nkt, selection.method = "vst", nfeatures = 2000)
K_ss_1_nkt[["mouse"]] <- "kalwrij"
K_ss_1_nkt[["timepoint"]] <- "SS"
K_ss_1_nkt[["number"]] <- "kalwrij_SS_1"

# Loading and pre-processing KalwRij SS 2 nk/t
K_ss_2_nkt <- Read10X(data.dir = "~/kal_2_ss/filtered_feature_bc_matrix/")
K_ss_2_nkt <- CreateSeuratObject(counts = K_ss_2_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=K_ss_2_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = K_ss_2_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = K_ss_2_nkt, slot = "counts"))
K_ss_2_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = K_ss_2_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
K_ss_2_nkt <- subset(x = K_ss_2_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_K_ss_2_nkt <- as.character(read.csv(file = "~/barcodes_kss2_nkt.csv", header = T, row.names = 1)[,1])
nk_K_ss_2_nkt <- gsub("_.", "", nk_K_ss_2_nkt)
K_ss_2_nkt <- subset(x=K_ss_2_nkt, cells = nk_K_ss_2_nkt)
K_ss_2_nkt <- NormalizeData(object = K_ss_2_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
K_ss_2_nkt <- FindVariableFeatures(object = K_ss_2_nkt, selection.method = "vst", nfeatures = 2000)
K_ss_2_nkt[["mouse"]] <- "kalwrij"
K_ss_2_nkt[["timepoint"]] <- "SS"
K_ss_2_nkt[["number"]] <- "kalwrij_SS_2"

# Loading and pre-processing KalwRij d21 3 nk/t
K_mm_3_nkt <- Read10X(data.dir = "~/kal_3_d21/filtered_feature_bc_matrix/")
K_mm_3_nkt <- CreateSeuratObject(counts = K_mm_3_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=K_mm_3_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = K_mm_3_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = K_mm_3_nkt, slot = "counts"))
K_mm_3_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = K_mm_3_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
K_mm_3_nkt <- subset(x = K_mm_3_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_K_mm_3_nkt <- as.character(read.csv(file = "~/barcodes_kmm3_nkt.csv", header = T, row.names = 1)[,1])
nk_K_mm_3_nkt <- gsub("_.", "", nk_K_mm_3_nkt)
K_mm_3_nkt <- subset(x=K_mm_3_nkt, cells = nk_K_mm_3_nkt)
K_mm_3_nkt <- NormalizeData(object = K_mm_3_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
K_mm_3_nkt <- FindVariableFeatures(object = K_mm_3_nkt, selection.method = "vst", nfeatures = 2000)
K_mm_3_nkt[["mouse"]] <- "kalwrij"
K_mm_3_nkt[["timepoint"]] <- "Diseased"
K_mm_3_nkt[["number"]] <- "kalwrij_d21_1"

# Loading and pre-processing KalwRij d21 4 nk/t
K_mm_4_nkt <- Read10X(data.dir = "~/kal_4_d21/filtered_feature_bc_matrix/")
K_mm_4_nkt <- CreateSeuratObject(counts = K_mm_4_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=K_mm_4_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = K_mm_4_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = K_mm_4_nkt, slot = "counts"))
K_mm_4_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = K_mm_4_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
K_mm_4_nkt <- subset(x = K_mm_4_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_K_mm_4_nkt <- as.character(read.csv(file = "~/barcodes_kmm4_nkt.csv", header = T, row.names = 1)[,1])
nk_K_mm_4_nkt <- gsub("_.", "", nk_K_mm_4_nkt)
K_mm_4_nkt <- subset(x=K_mm_4_nkt, cells = nk_K_mm_4_nkt)
K_mm_4_nkt <- NormalizeData(object = K_mm_4_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
K_mm_4_nkt <- FindVariableFeatures(object = K_mm_4_nkt, selection.method = "vst", nfeatures = 2000)
K_mm_4_nkt[["mouse"]] <- "kalwrij"
K_mm_4_nkt[["timepoint"]] <- "Diseased"
K_mm_4_nkt[["number"]] <- "kalwrij_d21_2"

# Loading and pre-processing BL6 SS 1 nk/t
B_ss_1_nkt <- Read10X(data.dir = "~/b6_1_ss/filtered_feature_bc_matrix/")
B_ss_1_nkt <- CreateSeuratObject(counts = B_ss_1_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=B_ss_1_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = B_ss_1_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = B_ss_1_nkt, slot = "counts"))
B_ss_1_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = B_ss_1_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
B_ss_1_nkt <- subset(x = B_ss_1_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_B_ss_1_nkt <- as.character(read.csv(file = "~/barcodes_bss1_nkt.csv", header = T, row.names = 1)[,1])
nk_B_ss_1_nkt <- gsub("_.", "", nk_B_ss_1_nkt)
B_ss_1_nkt <- subset(x=B_ss_1_nkt, cells = nk_B_ss_1_nkt)
B_ss_1_nkt <- NormalizeData(object = B_ss_1_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
B_ss_1_nkt <- FindVariableFeatures(object = B_ss_1_nkt, selection.method = "vst", nfeatures = 2000)
B_ss_1_nkt[["mouse"]] <- "BL6"
B_ss_1_nkt[["timepoint"]] <- "SS"
B_ss_1_nkt[["number"]] <- "BL6_SS_1"

# Loading and pre-processing BL6 SS 2 nk/t
B_ss_2_nkt <- Read10X(data.dir = "~/b6_2_ss/filtered_feature_bc_matrix/")
B_ss_2_nkt <- CreateSeuratObject(counts = B_ss_2_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=B_ss_2_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = B_ss_2_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = B_ss_2_nkt, slot = "counts"))
B_ss_2_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = B_ss_2_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
B_ss_2_nkt <- subset(x = B_ss_2_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_B_ss_2_nkt <- as.character(read.csv(file = "~/barcodes_bss2_nkt.csv", header = T, row.names = 1)[,1])
nk_B_ss_2_nkt <- gsub("_.", "", nk_B_ss_2_nkt)
B_ss_2_nkt <- subset(x=B_ss_2_nkt, cells = nk_B_ss_2_nkt)
B_ss_2_nkt <- NormalizeData(object = B_ss_2_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
B_ss_2_nkt <- FindVariableFeatures(object = B_ss_2_nkt, selection.method = "vst", nfeatures = 2000)
B_ss_2_nkt[["mouse"]] <- "BL6"
B_ss_2_nkt[["timepoint"]] <- "SS"
B_ss_2_nkt[["number"]] <- "BL6_SS_2"

# Loading and pre-processing BL6 d22 2 nk/t
B_mm_2_nkt <- Read10X(data.dir = "~/b6_2_ss/filtered_feature_bc_matrix/")
B_mm_2_nkt <- CreateSeuratObject(counts = B_mm_2_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=B_mm_2_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = B_mm_2_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = B_mm_2_nkt, slot = "counts"))
B_mm_2_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = B_mm_2_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
B_mm_2_nkt <- subset(x = B_mm_2_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_B_mm_2_nkt <- as.character(read.csv(file = "~/barcodes_bmm2_nkt.csv", header = T, row.names = 1)[,1])
nk_B_mm_2_nkt <- gsub("_.", "", nk_B_mm_2_nkt)
B_mm_2_nkt <- subset(x=B_mm_2_nkt, cells = nk_B_mm_2_nkt)
B_mm_2_nkt <- NormalizeData(object = B_mm_2_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
B_mm_2_nkt <- FindVariableFeatures(object = B_mm_2_nkt, selection.method = "vst", nfeatures = 2000)
B_mm_2_nkt[["mouse"]] <- "BL6"
B_mm_2_nkt[["timepoint"]] <- "Diseased"
B_mm_2_nkt[["number"]] <- "BL6_d22_2"

#Loading and pre-processing BL6 D22 3 nk/t
B_mm_3_nkt <- Read10X(data.dir = "~/b6_2_ss/filtered_feature_bc_matrix/")
B_mm_3_nkt <- CreateSeuratObject(counts = B_mm_3_nkt, min.cells = 3, min.features = 200, project = "mice")
mito.features_object2 <- grep(pattern = "^mt-", x=rownames(x=B_mm_3_nkt), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = B_mm_3_nkt, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = B_mm_3_nkt, slot = "counts"))
B_mm_3_nkt[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = B_mm_3_nkt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
B_mm_3_nkt <- subset(x = B_mm_3_nkt, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & nCount_RNA >200 & nCount_RNA <50000 & percent.mito <0.1)
nk_B_mm_3_nkt <- as.character(read.csv(file = "~/barcodes_bmm3_nkt.csv", header = T, row.names = 1)[,1])
nk_B_mm_3_nkt <- gsub("_.", "", nk_B_mm_3_nkt)
B_mm_3_nkt <- subset(x=B_mm_3_nkt, cells = nk_B_mm_3_nkt)
B_mm_3_nkt <- NormalizeData(object = B_mm_3_nkt, normalization.method = "LogNormalize", scale.factor = 1e4)
B_mm_3_nkt <- FindVariableFeatures(object = B_mm_3_nkt, selection.method = "vst", nfeatures = 2000)
B_mm_3_nkt[["mouse"]] <- "BL6"
B_mm_3_nkt[["timepoint"]] <- "Diseased"
B_mm_3_nkt[["number"]] <- "BL6_d22_3"

## Create Kalwrij Steady state nkt Object
# Identification of integration anchors
Kssnkt.list <- c(K_ss_1_nkt, K_ss_2_nkt)
Kssnkt.anchors <- FindIntegrationAnchors(object.list = Kssnkt.list, dims = 1:30)
Kssnkt.together <- IntegrateData(anchorset = Kssnkt.anchors, dims = 1:30)

## Create Kalwrij MM nkt Object
# Identification of integration anchors
Kmmnkt.list <- c(K_mm_1_nkt, K_mm_2_nkt)
Kmmnkt.anchors <- FindIntegrationAnchors(object.list = Kmmnkt.list, dims = 1:30)
Kmmnkt.together <- IntegrateData(anchorset = Kmmnkt.anchors, dims = 1:30)

## Create Kalwrij full nkt Object
# Identification of integration anchors
Knkt.list <- c(Kssnkt.together, Kmmnkt.together)
Knkt.anchors <- FindIntegrationAnchors(object.list = Knkt.list, dims = 1:30)
Knkt.together <- IntegrateData(anchorset = Knkt.anchors, dims = 1:30)

## Create B6 Steady state nkt Object
# Identification of integration anchors
Bssnkt.list <- c(B_ss_1_nkt, B_ss_2_nkt)
Bssnkt.anchors <- FindIntegrationAnchors(object.list = Bssnkt.list, dims = 1:30)
Bssnkt.together <- IntegrateData(anchorset = Bssnkt.anchors, dims = 1:30)

## Create B6 MM nkt Object
# Identification of integration anchors
Bmmnkt.list <- c(B_mm_2_nkt, B_mm_3_nkt)
Bmmnkt.anchors <- FindIntegrationAnchors(object.list = Bmmnkt.list, dims = 1:30)
Bmmnkt.together <- IntegrateData(anchorset = Bmmnkt.anchors, dims = 1:30)

## Create B6 full nkt Object
# Identification of integration anchors
Bnkt.list <- c(Bssnkt.together, Bmmnkt.together)
Bnkt.anchors <- FindIntegrationAnchors(object.list = Bnkt.list, dims = 1:30)
Bnkt.together <- IntegrateData(anchorset = Bnkt.anchors, dims = 1:30)

## Create full nkt Object
# Identification of integration anchors
nkt.list <- c(Knkt.together, Bnkt.together)
nkt.anchors <- FindIntegrationAnchors(object.list = Bnkt.list, dims = 1:30)
nkt.together <- IntegrateData(anchorset = Bnkt.anchors, dims = 1:30)

# Post-processing merged data
DefaultAssay(object = nkt.together) <- "integrated"
nkt.together <- ScaleData(object = nktl.together, verbose = F)
nkt.together <- RunPCA(object = nkt.together, verbose = F)
ElbowPlot(object = nkt.together, ndims = 40)
nkt.together <- RunUMAP(object = nkt.together, reduction = "pca", dims = 1:20)
nkt.together <- FindNeighbors(object = nkt.together, dims = 1:20)
nkt.together <- FindClusters(object = nkt.together, resolution = 0.3)
