#######IMPORTING WU DATASET, DATA QUALITY CONTROL, CLUSTERING AND SINGLE-CELL ANNOTATION WORKFLOW

####Loading necessary libraries.
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratData)
library(ggplot2)

###Loading single-cell data collected by Wu et al (2021). 
Wu.data<-Read10X("~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/filtered_gene_bc_matrices/Wu_etal_2021_BRCA_scRNASeq",
                 gene.column = 1)
Wu<-CreateSeuratObject(counts=Wu.data, project= "Wu3k",min.cells=3,min.features=200, meta.data = "metadata.csv")
Idents(Wu)=1
#Inspecting data properties.
Wu[["percent.mito"]]<-PercentageFeatureSet(Wu,pattern = "^MT")
VlnPlot(Wu,features = c("nCount_RNA","nFeature_RNA","percent.mito"),ncol = 3,alpha=0.25)
plot1<-FeatureScatter(Wu,feature1 = "nCount_RNA","percent.mito")
plot2<-FeatureScatter(Wu,feature1 = "nCount_RNA","nFeature_RNA")
plot1 + plot2
#Filtering out poor data, indicative of cell stress or artifacts from cell-sorting process.
Wu<-subset(Wu,subset=nFeature_RNA>200 & nFeature_RNA<4000 & percent.mito <5)
####Clustering workflow. 
#Normalizing data to account for non-biological variance.
Wu<-NormalizeData(Wu)
#Identification of highly variable features (unique genes).
Wu<-FindVariableFeatures(Wu,selection.method = "vst",nfeatures = 2000)
#Scaling data 
all.genes<-rownames(Wu)
Wu<-ScaleData(Wu,features=all.genes)
#Performing principal component analysis (PCA),reducing dimensions to ease data wrangling and handling.
Wu<-RunPCA(Wu,features=VariableFeatures(object=Wu))
#Inspecting dimension cut-off point. Choosing dimension cut-off point that is best trade-off
#for retaining most information, whilst making data wrangling and handling more feasible.
ElbowPlot(Wu)
#Choice of dimension is utilized for clustering data and reducing dimensions via 
#UMAP (Uniform Manifold Approximation and Projection). Resolution utilized for exemplar 
#run for cluster annotation. Deemed suitable for discerning cluster identity for the duration of project. 
Wu<-FindNeighbors(Wu,dims=1:15)
Wu<-FindClusters(Wu,resolution = 0.3)
Wu<-RunUMAP(Wu,dims=1:15)
#Producing DimPlot, displaying all 20 clusters identified.
DimPlot(Wu,reduction="umap")

####Single-cell annotation workflow.

###Inspecting markers for all clusters identified and producing gene heatmap.
Wu_markers<-FindAllMarkers(Wu,only.pos = TRUE)
#Identifying top 5 gene markers for exemplar gene heat map.
top5<-Wu_markers %>% 
  group_by(cluster) %>% 
  filter(avg_log2FC>1) %>% 
  slice_head(n=5) %>% 
  ungroup()
#Producing exemplar heat map.
DoHeatmap(Wu,features = top5$gene) + NoLegend() + theme(axis.text.y =element_text(size = 5)) +
  theme(panel.spacing.x = unit(0.5, "cm")) 

###Cluster idetntity interrogation
#Clusters were interrogated using violin plots, feature plots and differential expression
#using the FindMarkers() command that can be used for studying gene expression in clusters 
#or between clusters. It should be noted these annotations may not reflect the real identity of 
#cell clusters, being only an exemplar workflow. 

##Cluster 10 highly expresses MS4A1 relative to other clusters, likely representing 'B-cells'.
Cluster10<-FindMarkers(Wu,ident.1 = 10)
head(Cluster10,n=5)
VlnPlot(Wu,features=c("MS4A1"))
FeaturePlot(Wu,features = c("MS4A1"))

##Cluster 4 highly expresses RAMP2,GNG11 and PLVAP relative to other clusters. 
#Points to possible endothelial cell identity. Further interrogation shows high expression
#of PECAM1,VWF and EMCN. Therefore, cluster 4 is potentially of 'endothelial' identity. 
Cluster_4<-FindMarkers(Wu,ident.1 = 4)
head(Cluster_4,n=5)
VlnPlot(Wu,features = c("RAMP2","GNG11","PLVAP","VWF","PECAM1","EMCN"))
FeaturePlot(Wu,features = c("RAMP2","GNG11","PLVAP","VWF","PECAM1","EMCN"))

##Cluster 2 highly expresses LUM, CTSK, DCN and SFRP2, exclusively compared to other clusters. 
#Along with high COL1A2 expression. May point to cluster being of 'fibroblast' identity.
Cluster_2<-FindMarkers(Wu,ident.1=2)
head(Cluster_2,n=5)
VlnPlot(Wu,features = c("DCN","LUM","CTSK","COL1A2"))
FeaturePlot(Wu,features = c("DCN","LUM","CTSK","COL1A2","SFRP2"))

##Cluster 6 highly expresses IGFBP7, MYL9, NDUFA4L2 and ACTA2, whilst being negative for 
#Fibroblast associated markers like LUM,DCN or SFRP2. Likely identity is of 'smooth muscle cell'. 
Cluster_6<-FindMarkers(Wu,ident.1 = 6)
head(Cluster_6,n=5)
Cluster_6vs2<-FindMarkers(Wu,ident.1 = 6,ident.2 = 2)
VlnPlot(Wu,features = c("IGFBP7","MYL9","NDUFA4L2","ACTA2","SFRP2","LUM","DCN"))
FeaturePlot(Wu,features = c("IGFBP7","MYL9","SPARC","NDUFA4L2","SPARCL1","SFRP2","LUM","DCN"))

##Cluster 18 shows high exclusive expression of KRT14 and C2orf40, whilst a low expression of 
#EPCAM. Thus it likely points to a 'myoepithelial breast cell' identity. 
Cluster18<-FindMarkers(Wu,ident.1 = 18)
head(Cluster18,n=5)
VlnPlot(Wu,features = c("KRT14","C2orf40","EPCAM"))
FeaturePlot(Wu,features = c("KRT14","C2orf40","EPCAM"))

##Cluster 3 and 11 highly expresses macrophage associated markers like MS4A6A and CD68. 
#The largest gene expression difference is noted in FOLR2 expression, 
#highly up-regulated in cluster 11. Thus cluster 3 is FOLR2 negative and cluster 11 is FOLR2 positive. 
Cluster_3<-FindMarkers(Wu,ident.1 = 3)
head(Cluster_3,n=5)
Cluster_11<-FindMarkers(Wu,ident.1 = 11)
head(Cluster_11,n=5)
VlnPlot(Wu,features = c("MS4A6A","CD68"))
FeaturePlot(Wu,features=c("MS4A6A","CD68"))
Cluster_11vs3<-FindMarkers(Wu,ident.1 = 11,ident.2 = 3) 
head(Cluster_11vs3,n=5) 
VlnPlot(Wu,features = c("FOLR2"))
FeaturePlot(Wu,features = c("FOLR2"))

##Clusters 5,14 and 17 are enriched for MZB1,DERL3,JCHAIN. Clusters 5 and 14 highly express immunoglobulin-related genes.
#This points to their plasma cell identity. Clusters 5 and 14 do not substantially differ from one another
#therefore both labelled as plasma cells. Cluster 17 displays particularly high expression of GZMB,
#LILRA4 and PLAC8. This may point to a potential 'plasmacytoid dendritic cell (pDC)' identity.
Cluster5<-FindMarkers(Wu,ident.1 = 5)
head(Cluster5,n=5)
Cluster14<-FindMarkers(Wu,ident.1=14)
head(Cluster14,n=5)
VlnPlot(Wu,features = c("MZB1","JCHAIN","IGHG1","IGHG4","DERL3"))
FeaturePlot(Wu,features = c("MZB1","JCHAIN","IGHG1","IGHG4","DERL3"))
Cluster_17<-FindMarkers(Wu,ident.1 = 17)
head(Cluster_17,n=5)
Cluster_17vs5and4<-FindMarkers(Wu,ident.1 = 17,ident.2 = c(5,4))
head(Cluster_17vs5and4,n=5)
VlnPlot(Wu,features=c("PLAC8","LILRA4","GZMB","IRF8"))
FeaturePlot(Wu,features=c("PLAC8","LILRA4","GZMB","IRF8"))

##Clusters 0,1,9,12 and 19 display higher CD3 expression compared to other clusters, pointing to potential T-cell clusters. 
#Nevertheless, discerning their phenotypes is complex. Thus, these annotations are only for discernment.
VlnPlot(Wu,features = c("CD3D","CD3E"))
FeaturePlot(Wu,features = c("CD3D","CD3E"))

#Cluster 0 shows a rich population of IL7R cells relative to other T-cell clusters, with low CD4 and CD8 expression. 
#Due to complexity of phenotype assigned as 'IL7R rich T-cells'. 
Cluster_0<-FindMarkers(Wu,ident.1 = 0)
head(Cluster_0,n=5)
Cluster_0vsall<-FindMarkers(Wu,ident.1 = 0,ident.2 = c(1,9,12,19))
head(Cluster_0vsall,n=5)
VlnPlot(Wu,features = c("IL7R","CD3D","CD8A","CD4"))
FeaturePlot(Wu,features = c("IL7R","CD3D","CD8A","CD4"))

#Cluster 1 displays the highest expression of NKG7, CCL4, CCL5 and GZMB, relative to other T-cell clusters. 
#Thus assigned as 'NKG7 rich T-cells'. 
Cluster_1<-FindMarkers(Wu,ident.1 = 1)
head(Cluster_1,n=5)
Cluster_1vsall<-FindMarkers(Wu,ident.1 = 1,ident.2=c(0,9,12,19))
head(Cluster_1vsall,n=5)
VlnPlot(Wu,features = c("NKG7","CCL4","CCL5","GZMB"))
FeaturePlot(Wu,features = c("NKG7","CCL4","CCL5","GZMB"))

#Cluster 9 shows high expression of CD3E and a plethora of immunoglobulin-genes relative to other T-cells
#Thus, assigned as 'IGKV1-5 positive, T-cell rich'. 
Cluster_9<-FindMarkers(Wu,ident.1 = 9)
head(Cluster_9,n=5)
Cluster_9vsall<-FindMarkers(Wu,ident.1 = 9,ident.2 = c(0,1,12,19))
head(Cluster_9vsall,n=5)
VlnPlot(Wu,features = c("IGKV1-5","IGKV1-9","IGKV3-11","IGKV3-20","IGLV3-25"))
FeaturePlot(Wu,features = c("IGKV1-5","IGKV1-9","IGKV3-11","IGKV3-20","IGLV3-25"))

#Cluster 12 shows high expression of STMN1,MKI67,TYMS,NUSAP and CKS1B, which are markers of proliferation. 
#Cluster 12 also displays high expression of CD3. Thus, assigned as 'cycling cells, T-cell rich'.
Cluster_12<-FindMarkers(Wu,ident.1 = 12)
head(Cluster_12,n=5)
Cluster_12vsall<-FindMarkers(Wu,ident.1 = 12,ident.2 = c(0,1,9,19))
head(Cluster_12vsall,n=5)
VlnPlot(Wu,features = c("STMN1","MKI67","TYMS","NUSAP1","CKS1B","CD3E"))
FeaturePlot(Wu,features = c("STMN1","MKI67","TYMS","NUSAP1","CKS1B","CD3E"))

#Cluster 19 possess high expression of NKG7 wit high S100A8 and S100A9 expression.
#Furthermore, GZM genes are highly up-regulated with some CD3 and 
#Killer cell receptor (KLR) positive cells. May point to potential 'NK positive, T cell positive' 
Cluster_19<-FindMarkers(Wu,ident.1 = 19)
head(Cluster_19,n=5)
Cluster_19vsall<-FindMarkers(Wu,ident.1 = 19,ident.2 = c(0,1,9,12))
head(Cluster_19vsall,n=5)
VlnPlot(Wu,features = c("GZMB","S100A9","S100A8","KLRC2","CD3D","KLRD1","NKG7"))
FeaturePlot(Wu,features = c("GZMB","S100A9","S100A8","KLRC2","CD3D","KLRD1","NKG7","ZMAT4"))

##Clusters 7,8,13,15 and 16 express EPCAM,KR substantially higher compared to other clusters
#pointing to potential epithelial cell identity.
VlnPlot(Wu,features = c("EPCAM"))
FeaturePlot(Wu,features = c("EPCAM"))

#Cluster 7, has shown to express CRYAB,SAA1 and KRT15, exclusively relative to other epithelial clusters
#Largest exclusive expression was for CRYAB and SAA1. Thus cluster assigned as 
#'CRYAB,SAA1 positive epithelial'. 
Cluster_7<-FindMarkers(Wu,ident.1=7)
head(Cluster_7,n=5)
Cluster_7vsall<-FindMarkers(Wu,ident.1 = 7,ident.2 = c(8,13,15,16))
Cluster_7vsall$GeneName<-rownames(Cluster_7vsall)
head(Cluster_7vsall,n=5)
FeaturePlot(Wu,features = c("KRT15","CRYAB","SAA1"))
VlnPlot(Wu,features = c("KRT15","CRYAB","SAA1"))

#Cluster 8, has shown to express AGR3,ANKRD30A and MLPH, exclusively relative to other epithelial clusters
#Largest expression was shown for AGR3, thus cluster assigned as 'AGR3 positive epithelial'. 
Cluster_8<-FindMarkers(Wu,ident.1=8)
head(Cluster_8,n=5)
Cluster_8vsall<-FindMarkers(Wu,ident.1 = 8,ident.2 = c(7,13,15,16))
head(Cluster_8vsall,n=5)
VlnPlot(Wu,features = c("AGR3","MLPH","ANKRD30A"))
FeaturePlot(Wu,features = c("AGR3","ANKRD30A","MLPH"))

#Cluster 13 expresses high levels of S100A8,S100A9 and KRT17, exclusively relative to other epithelial clusters. 
#Thus assigned as 'S100A8/9,KRT17 positive epithelial'.
Cluster_13<-FindMarkers(Wu,ident.1 = 13)
head(Cluster_13,n=5)
Cluster_13vsall<-FindMarkers(Wu,ident.1 = 13,ident.2 = c(7,8,15,16))
head(Cluster_13vsall,n=5)
VlnPlot(Wu,features = c("S100A8","S100A9","KRT17"))
FeaturePlot(Wu,features = c("S100A8","S100A9","KRT17"))

#Cluster 15 expresses high levels of FDCSP, exclusively relative to other epithelial clusters.
#Thus assigned as 'FDCSP positive epithelial' 
Cluster_15<-FindMarkers(Wu,ident.1 = 15)
head(Cluster_15,n=5)
Cluster_15vsall<-FindMarkers(Wu,ident.1 = 15,ident.2 = c(7,8,13,16))
head(Cluster_15vsall,n=5)
VlnPlot(Wu,features = c("FDCSP"))
FeaturePlot(Wu,features=c("FDCSP"))

#Cluster 16 expresses high levels of SCGB1B2P,MUCL1 and SCGB2B2 relative to other epithelial cluster
#Highest expression is of SCGB1B2P, thus assigned as 'SCGB1B2P positive epithelial'. 
Cluster_16<-FindMarkers(Wu,ident.1=16)
head(Cluster_16,n=5)
Cluster_16vsall<-FindMarkers(Wu,ident.1 = 16,ident.2 = c(7,8,13,15))
head(Cluster_16vsall,n=5)
VlnPlot(Wu,features=c("SCGB1B2P","MUCL1","SCGB2B2"))
FeaturePlot(Wu,features=c("SCGB1B2P","MUCL1","SCGB2B2"))



###Assigning cluster names accordingly. 
#Naming clusters according to cell type identified
New.cluster.ids1<-c("IL7R rich T-cells","NKG7 rich T-cells","Fibroblasts","FOLR2 negative macrophages",
                    "Endothelial cells","Plasma cells","Smooth muscle cells","CRYAB,SAA1 positive epithelial cells",
                    "AGR3 positive epithelial cells","IGKV1-5 positive, T-cell rich","B cells","FOLR2 positive macrophages",
                    "Cycling cells, T-cell rich","S100A8/9,KRT17 positive epithelial cells","Plasma cells",
                    "FDCSP positive epithelial cells","SCGB1B2P positive epithelial cells",
                    "pDC","Myoepithelial cells","NK positive, T-cell positive")
names(New.cluster.ids1)<-levels(Wu)
Wu<-RenameIdents(Wu,New.cluster.ids1)
DimPlot(Wu,reduction="umap",label=TRUE,pt.size=0.05,label.box = TRUE,label.size = 2.5,repel = TRUE) + NoLegend()




