####RCTD DECONVOLUTION OF BREAST TISSUE SAMPLE USIGN SELF-MADE SINGLE-CELL ANNOTATIONS.
###LOADING NECESSARY LIBRARIES
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(spacexr)
library(Matrix)
library(doParallel)
###Importing breast image, spatial transcriptomics data and subjugating it to single-cell clustering 
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
##Carrying out typical workflow of pre-processing spatial transcriptomics data (see:https://satijalab.org/seurat/articles/spatial_vignette) 
Breast<-NormalizeData(Breast,verbose=FALSE,assay="Spatial")
Breast<-FindVariableFeatures(Breast,assay = "Spatial")
all.genes<-rownames(Breast)
Breast<-ScaleData(Breast,features=all.genes)
Breast<-RunPCA(Breast,verbose=FALSE,assay="Spatial")
ElbowPlot(Breast)
Breast <- FindNeighbors(Breast, reduction = "pca", dims = 1:13)
Breast <- FindClusters(Breast, verbose = FALSE)
Breast <- RunUMAP(Breast, reduction = "pca", dims = 1:13)
###MAKING SINGLE-CELL REFERENCES BASED ON ANNOTATIONS. REQUIRES TWO INPUTS: 'COUNTS' and 'CELL_TYPES'
##Counts is a matrix representing the Digital Gene Expression
counts<-Wu.data
##Cell_types are the named factor of cell type for each cell annotated.
#Obtaining cell types from the self-annotated data
meta_data<-Wu@meta.data
meta_data<- tibble::rownames_to_column(meta_data, "X") 
cell_types<-meta_data$renamedCellTypes
names(cell_types)<-meta_data$X
cell_types<-as.factor(cell_types)
#Constructing a single-cell reference
reference<-Reference(counts,cell_types)
print(dim(reference@counts))
table(reference@cell_types)
###Loading Spatial Transcriptomics data as a 'SpatialRNA' object. 
##Creating SpatialRNA object
SpatialRNA<-read.VisiumSpatialRNA("~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/Visualization")
###Creating and running RCTD in full 'doubletmode' for single-cell annotations. 
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
##Cell contribution of spots mapped onto the original breast tissue biological image. 
#Making a dataframe with cell weights
newbarcode<-norm_weights@Dimnames[[1]]
nUMIs<-norm_weights@Dimnames[[1]]

newdata <- data.frame(Bcells_weight = norm_weights[,'B cells'],
                      Fibroblast_weight=norm_weights[,'Fibroblasts'],
                      NKG7_Tcell_weight=norm_weights[,'NKG7 rich T-cells'],
                      IL7R_Tcell_weight=norm_weights[,'IL7R rich T-cells'],
                      Cycling_cells_Tcell_weight=norm_weights[,'Cycling cells, T-cell rich'],
                      Endothelial_weight=norm_weights[,'Endothelial cells'],
                      pDC_weight=norm_weights[,'pDC'],FOLR2_negative_macrophage_weight=norm_weights[,'FOLR2 negative macrophages'],
                      FOLR2_positive_macrophage_weight=norm_weights[,'FOLR2 positive macrophages'],
                      Smoothmuscle_weight=norm_weights[,'Smooth muscle cells'],
                      CRYAB_SAA1_positive_epithelial_weight=norm_weights[,'CRYAB,SAA1 positive epithelial cells'],
                      FDCSP_positive_epithelial_cell_weight=norm_weights[,'FDCSP positive epithelial cells'],
                      AGR3_positive_epithelial_cell_weight=norm_weights[,'AGR3 positive epithelial cells'],
                      Plasma_cells_weight=norm_weights[,'Plasma cells'],
                      NKcell_Tcell_weight=norm_weights[,'NK cell positive, T-cell positive'],
                      Myoepithelial_weight=norm_weights[,'Myoepithelial cells'],
                      S100A8_KRT17_epithelial_cell_weight=norm_weights[,'S100A8,KRT17 positive epithelial cells'],
                      IGKV15_Tcell_weight=norm_weights[,'IGKV1-5 positive, T-cell rich'],
                      SCGB1B2P_epithelial_weight=norm_weights[,'SCGB1B2P positive epithelial cells'])

#Add the new metadata to breast image
Breast_integrated <- AddMetaData(Breast, metadata = newdata)

###RCTD plots
#B-cells weight

Breast_integrated$'Bcell weight' = newdata$Bcells_weight
Bcellplot<-SpatialFeaturePlot(Breast_integrated,"Bcell weight",
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15),
  )
Bcellplot

#Fibroblast weight
Breast_integrated$'Fibroblast weight'=newdata$Fibroblast_weight
Fibroblastplot<-SpatialFeaturePlot(Breast_integrated,"Fibroblast weight",
                            image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Fibroblastplot

#Myoepithelial weight
Breast_integrated$'Myoepithelial weight'=newdata$Myoepithelial_weight
Myoepithelialplot<-SpatialFeaturePlot(Breast_integrated,"Myoepithelial weight",
                                   image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
Myoepithelialplot

#FOLPR2 positive macrophage weight
Breast_integrated$'FOLPR2 positive macrophage weight'=newdata$FOLR2_positive_macrophage_weight
FOLPR2_pos<-SpatialFeaturePlot(Breast_integrated,'FOLPR2 positive macrophage weight',
                                    image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
FOLPR2_pos

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












































































