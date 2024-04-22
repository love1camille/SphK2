library(tidyverse)
library(Seurat)
library(ggplot2)
library(clustree)
library (cowplot)
library (dplyr)
library(patchwork)

setwd("GSE182365")
CH <- Read10X(data.dir ="./CH")
CH <- CreateSeuratObject(counts = CH, project = "CH",min.cells = 3, min.features = 200)

HH <- Read10X(data.dir ="./HH")
HH <- CreateSeuratObject(counts = HH, project = "HH",min.cells = 3, min.features = 200)

CH[["percent.mt"]] <- PercentageFeatureSet(CH,pattern = "^MT-")
HH[["percent.mt"]] <- PercentageFeatureSet(HH,pattern = "^MT-")
preQC_CH <- VlnPlot(CH, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    group.by = "orig.ident", 
                    pt.size = 0)
preQC_HH <- VlnPlot(HH, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, 
                    group.by = "orig.ident", 
                    pt.size = 0)

CH <- subset(CH, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 15)
HH <- subset(HH, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 15)
postQC_CH <- VlnPlot(CH, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     ncol = 3, 
                     group.by = "orig.ident", 
                     pt.size = 0)
postQC_HH <- VlnPlot(HH, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                     ncol = 3, 
                     group.by = "orig.ident", 
                     pt.size = 0)
postQC_CH+postQC_HH

HH <- NormalizeData(HH)
HH <- FindVariableFeatures(HH, nfeatures = 2000)
CH <- NormalizeData(CH)
CH <- FindVariableFeatures(CH, nfeatures = 2000)


sampleList <- list(CH, HH)
H <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
H <- IntegrateData(anchorset = H, dims = 1:50)
save(H, file = "H.RData")

H <- ScaleData(H, verbose = FALSE)
H <- RunPCA(H, npcs = 50, verbose = FALSE)
H <- FindNeighbors(H, reduction = "pca", dims = 1:50)
H <- FindClusters(H, 
                  resolution = 0.05)
H <- RunUMAP(H, reduction = "pca", dims = 1:50)


DimPlot(H,label = T,split.by = "orig.ident",ncol = 3)

save(H, file = "H.RData")

DefaultAssay(H) <- "RNA"
all.markers  <- FindAllMarkers(H, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold =0.75)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.2, ]
write.csv(significant.markers, file = "significant.markers.csv")

markers <- c("Sphk2")
DotPlot(H,features = markers)+coord_flip()

FeaturePlot(H,features = c("Sphk2"), split.by = "orig.ident")

alldata <- ScaleData(H, 
                     features = markers, 
                     assay = "RNA")
DoHeatmap(alldata, 
          features = markers,
          group.by = "orig.ident",
          assay = "RNA")

FeaturePlot(H, features = markers, split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(H, features = "Klf4", split.by = "orig.ident", group.by = "seurat_clusters", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)