#######IMPORTING WU DATASET AND BREAST TISSUE SAMPLE FROM 10X GENOMICS AND CARRYING OUT EXEMPLAR RCTD OUTPUT.
#####Importing breast tissue sample from 10X genomics (Link to repository:https://www.10xgenomics.com/datasets/fresh-frozen-visium-on-cytassist-human-breast-cancer-probe-based-whole-transcriptome-profiling-2-standard)
setwd("~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/Visualization")
Data.dir<-"~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/Visualization"
list.files(Data.dir)
image.dir<-"~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/Visualization"
Breast<-Load10X_Spatial(data.dir=Data.dir,
                        filename="filtered_feature_bc_matrix.h5",
                        assay = "Spatial",
                        slice="slice1",
                        filter.matrix = TRUE,
                        to.upper = FALSE,
                        image=NULL)
###Carrying out typical workflow of pre-processing spatial transcriptomics data (see:https://satijalab.org/seurat/articles/spatial_vignette) 
Breast<-NormalizeData(Breast,verbose=FALSE,assay="Spatial")
Breast<-FindVariableFeatures(Breast,assay = "Spatial")
all.genes<-rownames(Breast)
Breast<-ScaleData(Breast,features=all.genes)
Breast<-RunPCA(Breast,verbose=FALSE,assay="Spatial")
ElbowPlot(Breast)
Breast <- FindNeighbors(Breast, reduction = "pca", dims = 1:13)
Breast <- FindClusters(Breast, verbose = FALSE)
Breast <- RunUMAP(Breast, reduction = "pca", dims = 1:13)
##Visualizing clusters
Bplot1<-DimPlot(Breast,reduction="umap",label = TRUE)
Bplot2<-SpatialDimPlot(Breast,label=TRUE,label.size=3)
Bplot1+Bplot2


####Importing single-cell data collected by Wu et al (2021)
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
###Carrying out typical workflow of pre-processing single-cell data before input into deconvolution tools. 
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
#run for cluster annotation. Deemed suitable for discerning cluster identity. 
Wu<-FindNeighbors(Wu,dims=1:15)
Wu<-FindClusters(Wu,resolution = 0.3)
Wu<-RunUMAP(Wu,dims=1:15)
#Producing DimPlot, displaying all 20 clusters identified.
DimPlot(Wu,reduction="umap")

####RCTD deconvolution of breast tissue sample using single-cell data annotated by Wu.
###Integrating breast tissue sample and Wu annotations 
anchors <- FindTransferAnchors(reference = Wu,
                               query = Breast, normalization.method = "LogNormalize")
predictions.assay <- TransferData(anchorset = anchors, refdata = Wu$seurat_clusters, prediction.assay = TRUE,
                                  weight.reduction = Breast[["pca"]], dims = 1:13)
Breast[["predictions"]] <- predictions.assay
DefaultAssay(Breast) <- "predictions"
SpatialFeaturePlot(Breast, features = c("0","1","2","3","4","5","6","7","8","9","10","11","12"), pt.size.factor = 1.6, ncol = 5, crop = TRUE)

####Running RCTD
#Installing the necessary libraries
library(spacexr)
library(Matrix)
library(doParallel)
### Making single-cell references. It requires two compulsory inputs: 'counts' and 'cell_types'.
##Counts is a matrix representing the Digital Gene Expression
counts<-Wu.data
##Cell_types are the named factor of cell type for each cell.
#Retrieving meta data of Wu scRNA data
meta_data<-read.csv(file.path(
  "~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/filtered_gene_bc_matrices/Wu_etal_2021_BRCA_scRNASeq","metadata.csv"))
#Obtaining cell type information 
cell_types<-meta_data$celltype_major
names(cell_types)<-meta_data$X
cell_types<-as.factor(cell_types)
cell_types
#Constructing a single-cell reference
reference<-Reference(counts,cell_types)
print(dim(reference@counts))
table(reference@cell_types)

###Loading Spatial Transcriptomics data as a 'SpatialRNA' object. 
##Creating SpatialRNA object
SpatialRNA<-read.VisiumSpatialRNA("~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/Visualization")

###Creating and running RCTD in full 'doubletmode' for 'major_cell_type' category within the Wu data. 
myRCTD<-create.RCTD(SpatialRNA,reference,max_cores=6)
myRCTD <- run.RCTD(myRCTD,doublet_mode = "full")
#Analyzing RCTD results
results<-myRCTD@results
#Normalizing cell type proportion to sum to 1
norm_weights=normalize_weights(results$weights)
#List of cell type names
cell_type_names<-myRCTD@cell_type_info$info[[2]]
spatialRNA<-myRCTD@spatialRNA
resultsdir<-'RCTD_plots'
dir.create(resultsdir)

###RCTD results
#Plots confident weights for each cell type as in full_mode
plot_weights(cell_type_names,spatialRNA,resultsdir,norm_weights)
#Plots weights for each cell type as in full_mode.
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
#Plots the number of confident pixels of each cell type in full_mode. 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

##Cell contribution of spots mapped onto the original breast tissue biological image. 
#Making a dataframe with cell weights
newbarcode<-norm_weights@Dimnames[[1]]
nUMIs<-norm_weights@Dimnames[[1]]
newdata <- data.frame(Bcells_weight = norm_weights[,'B-cells'],CAFs_weight = norm_weights[, 'CAFs'],
                      Cancerepithelial_weight = norm_weights[, 'Cancer Epithelial'], 
                      Endothelial_weight = norm_weights[, 'Endothelial'],  
                      Myeloid_weight = norm_weights[, 'Myeloid'],  
                      NormalEpithelial_weight = norm_weights[, 'Normal Epithelial'], 
                      Plasmablast_weight = norm_weights[, 'Plasmablasts'], 
                      PVL_weight = norm_weights[, 'PVL'], 
                      Tcells_weight = norm_weights[, 'T-cells'])
newdata_longer<-pivot_longer(newdata,cols=c(Bcells_weight,CAFs_weight,Cancerepithelial_weight,Endothelial_weight,
                                            Myeloid_weight,NormalEpithelial_weight,Plasmablast_weight,PVL_weight,Tcells_weight),
                             names_to = "cell_type",
                             values_to = "weight")
#Add the new metadata to breast image
Breast_integrated <- AddMetaData(Breast, metadata = newdata)

#B-cells weight

Breast_integrated$'Bcell weight' = newdata$Bcells_weight
Bcellplot<-SpatialFeaturePlot(Breast_integrated,"Bcell weight",
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15),
  )
Bcellplot

#CAFs weight
Breast_integrated$'CAFs weight'=newdata$CAFs_weight
CAFplot<-SpatialFeaturePlot(Breast_integrated,"CAFs weight",
                            image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
CAFplot

#Cancer epithelial weight
Breast_integrated$'Cancerepithelial weight'=newdata$Cancerepithelial_weight
Epithelialplot<-SpatialFeaturePlot(Breast_integrated,"Cancerepithelial weight",
                                   image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Epithelialplot

#Endothelial weight
Breast_integrated$'Endothelial weight'=newdata$Endothelial_weight
Endothelialplot<-SpatialFeaturePlot(Breast_integrated,'Endothelial weight',
                                    image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Endothelialplot

#Myeloid weight
Breast_integrated$'Myeloid weight'=newdata$Myeloid_weight
Myeloidplot<-SpatialFeaturePlot(Breast_integrated,'Myeloid weight',
                                image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Myeloidplot

#Normal epithelial weight
Breast_integrated$'NormalEpithelial weight'=newdata$NormalEpithelial_weight
Normalepithelialplot<-SpatialFeaturePlot(Breast_integrated,'NormalEpithelial weight',
                                         image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Normalepithelialplot

#Plasmablasts weight
Breast_integrated$'Plasmablast weight'=newdata$Plasmablast_weight
Plasmablastplot<-SpatialFeaturePlot(Breast_integrated,'Plasmablast weight',image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Plasmablastplot

#PVL_weight
Breast_integrated$'PVL weight'=newdata$PVL_weight
PVLplot<-SpatialFeaturePlot(Breast_integrated,'PVL weight',
                            image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
PVLplot

#Tcells_weight
Breast_integrated$'Tcell weight'=newdata$Tcells_weight
Tcellplot<-SpatialFeaturePlot(Breast_integrated,'Tcell weight',
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Tcellplot

#Combined plot of cell types weights
Image<-SpatialFeaturePlot(Breast_integrated,features=c('Bcell weight','CAFs weight',
                                                       'Cancerepithelial weight',
                                                       'Endothelial weight',
                                                       'Myeloid weight',
                                                       'NormalEpithelial weight',
                                                       'Plasmablast weight',
                                                       'PVL weight',
                                                       'Tcell weight'),ncol=5,
                          image.alpha = 0.5,alpha=c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=7.5),
        plot.title = element_text(size=14,face="bold",vjust=1,hjust=0.5))
Image













































































