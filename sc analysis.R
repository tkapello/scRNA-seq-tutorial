# Remove unneeded objects from memory
rm(list=ls())

# Load packages
library(Matrix)
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Set directory
counts<-readRDS('raw_counts.rds')
metadata<-read.csv('metadata.csv')

# Load data
seurat<-CreateSeuratObject(counts=as.matrix(counts), min.cells=5, min.features=80)

# Decide on inflection point
df<-data.frame(row.names=rownames(seurat@meta.data),
               cell=rownames(seurat@meta.data),
               Gene=seurat$nFeature_RNA,
               UMI=seurat$nCount_RNA)

df<-df[order(df$Gene, decreasing=T),]
df$index<-1:nrow(df)

ggplot(df, aes(index, Gene))+
  geom_point(size=0.1, show.legend=F)+
  scale_y_log10()+
  scale_x_log10()+
  geom_hline(aes(yintercept=200))+
  geom_hline(aes(yintercept=100), color='red')+
  geom_hline(aes(yintercept=50), color='blue')+
  ylab('UMI numbers')+
  xlab('cells sorted by gene count')+
  theme_bw()

ggplot(df, aes(Gene, UMI))+
  geom_point(size=0.1, show.legend=F)+
  scale_y_log10()+
  scale_x_log10()+
  geom_abline(slope=1, color='red')+
  geom_vline(aes(xintercept=100), color='orange')+
  geom_vline(aes(xintercept=200), color='brown')+
  ylab('UMI numbers')+
  xlab('Gene numbers')+
  theme_bw()

# Filter on number of genes
quantile(seurat$nFeature_RNA, probs=c(0.25, 0.5, 0.75, 0.95, 0.99))
seurat<-subset(seurat, subset=nFeature_RNA>=200)

# Filter on mitochondrial gene content
seurat[["percent.mt"]]<-PercentageFeatureSet(object=seurat, pattern="^mt-")
quantile(seurat$percent.mt, probs=(seq(0, 1, by=0.25)))

seurat<-subset(seurat, subset=percent.mt<=10)

# QC analysis
Idents(seurat)<-'orig.ident'

par(mfrow=c(1, 2))
FeatureScatter(object=seurat, feature1="nCount_RNA", feature2="percent.mt")
FeatureScatter(object=seurat, feature1="nFeature_RNA", feature2="nCount_RNA")

VlnPlot(object=seurat, features=c("nFeature_RNA", "nCount_RNA","percent.mt"), log=F, ncol=3)

# Normalisation 
seurat<-NormalizeData(object=seurat, normalization.method="LogNormalize", scale.factor=10000)

# Find variable genes
seurat<-FindVariableFeatures(object=seurat, selection.method="vst", nfeatures=2000)

# Scale data
seurat<-ScaleData(object=seurat)

# Perform linear reduction    
seurat<-RunPCA(object=seurat, features=seurat@assays$RNA@var.features, verbose=F, npcs=50)
print(x=seurat[["pca"]], dims=1:15, nfeatures=5)

VizDimLoadings(object=seurat, dims=1:15, reduction="pca")

DimPlot(object=seurat, reduction="pca")

DimHeatmap(object=seurat, dims=1:10, cells=500, balanced=TRUE)

# Decide on which PCs to keep
ElbowPlot(object=seurat, ndims=50)

# Find k-nearest neighbours
seurat<-FindNeighbors(object=seurat, reduction='pca', dims=1:20)

# Run UMAP
seurat<-RunUMAP(object=seurat, reduction="pca", dims=1:20, verbose=F)
DimPlot(object=seurat, reduction='umap', label=T)

# Find clusters
seurat<-FindClusters(object=seurat, resolution=0.8, group.singletons=T)

# Look for potential batches
DimPlot(seurat, group.by='orig.ident') + scale_color_jco()

# Calculate cluster markers
deg<-FindAllMarkers(seurat, min.pct=0.20, logfc.threshold=0.25, test.use='MAST', only.pos=T)
deg<-deg[deg$p_val_adj<=0.05,]
deg<-deg %>% arrange(cluster, -avg_log2FC)
deg5<-deg %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)

DotPlot(object=seurat, features=unique(deg5$gene), cols='RdBu', dot.scale=6, dot.min=0.1) + RotatedAxis()
