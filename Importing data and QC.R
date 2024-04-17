####IMPORTING WU DATASET AND CARRYING OUT DATA QUALITY CONTROL 
#Loading necessariy libraries
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratData)
library(ggplot2)
#Importing raw Wu dataset 
Wu.data<-Read10X("~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/filtered_gene_bc_matrices/Wu_etal_2021_BRCA_scRNASeq",
                 gene.column = 1)
#Creating a Seurat object and inspecting data
Wu<-CreateSeuratObject(counts=Wu.data, project= "Wu3k",min.cells=3,min.features=200, meta.data = "metadata.csv")
Idents(Wu)=1
Wu[["percent.mito"]] <- PercentageFeatureSet(Wu,pattern = "^MT-")
VlnPlot(Wu,features = c("nCount_RNA","nFeature_RNA","percent.mito"),ncol = 3,alpha=0.25)
plot1<-FeatureScatter(Wu,feature1 = "nCount_RNA","percent.mito")
plot2<-FeatureScatter(Wu,feature1 = "nCount_RNA","nFeature_RNA")
plot1 + plot2
#Filtering out poor data
Wu<-subset(Wu.1,subset=nFeature_RNA>200 & nFeature_RNA<4000 & percent.mito <5)
