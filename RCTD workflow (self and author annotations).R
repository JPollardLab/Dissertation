####RCTD DECONVOLUTION OF BREAST TISSUE SAMPLE USING SELF-MADE SINGLE-CELL ANNOTATIONS.
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
                      CAF_weight=norm_weights[,'CAFs'],
                      NKG7_Tcell_weight=norm_weights[,'NKG7 rich T-cells'],
                      IL7R_Tcell_weight=norm_weights[,'IL7R rich T-cells'],
                      Cycling_cells_Tcell_weight=norm_weights[,'Cycling cells, T-cell rich'],
                      Endothelial_weight=norm_weights[,'Endothelial cells'],
                      pDC_weight=norm_weights[,'pDC'],FOLR2_negative_macrophage_weight=norm_weights[,'FOLR2 negative macrophages'],
                      FOLR2_positive_macrophage_weight=norm_weights[,'FOLR2 positive macrophages'],
                      Smoothmuscle_weight=norm_weights[,'Smooth muscle cells'],
                      CRYAB_SAA1_positive_epithelial_weight=norm_weights[,'CRYAB,SAA1 positive epithelial cells'],
                      Cancer_epithelial_cell_weight=norm_weights[,'Cancer epithelial cells'],
                      Plasmablasts_weight=norm_weights[,'Plasmablasts'],
                      NKcell_Tcell_weight=norm_weights[,'NK cell positive, T-cell positive'],
                      Myoepithelial_weight=norm_weights[,'Myoepithelial cells'],
                      S100A8_epithelial_cell_weight=norm_weights[,'S100A8 positive epithelial cells'],
                      IGKV15_Tcell_weight=norm_weights[,'IGKV1-5 positive, T-cell rich'],
                      SCGB1B2P_epithelial_weight=norm_weights[,'SCGB1B2P positive epithelial cells'])

#Add the new metadata to breast image
Breast_integrated <- AddMetaData(Breast, metadata = newdata)

###RCTD plots
#B-cells weight

Breast_integrated$'Bcell weight self-annotation' = newdata$Bcells_weight
Bcellplot<-SpatialFeaturePlot(Breast_integrated,"Bcell weight self-annotation",
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Bcellplot

#CAF weight
Breast_integrated$'CAF weight self-annotation'=newdata$CAF_weight
CAF_plot<-SpatialFeaturePlot(Breast_integrated,"CAF weight self-annotation",
                                   image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=15))
CAF_plot


#NKG7 rich T-cell weight
Breast_integrated$'NKG7 rich T-cell weight self-annotation'=newdata$NKG7_Tcell_weight
NKG7_plot<-SpatialFeaturePlot(Breast_integrated,"NKG7 rich T-cell weight self-annotation",
                            image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
NKG7_plot

#IL7R rich T-cells
Breast_integrated$'IL7R rich T-cell weight self-annotation'=newdata$IL7R_Tcell_weight
IL7R_plots<-SpatialFeaturePlot(Breast_integrated,"IL7R rich T-cell weight self-annotation",
                                      image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
IL7R_plots

#Cycling cells, T-cell rich
Breast_integrated$'Cycling cells, T-cell rich weight self-annotation'=newdata$Cycling_cells_Tcell_weight
Cycling_plot<-SpatialFeaturePlot(Breast_integrated,"Cycling cells, T-cell rich weight self-annotation",
                               image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=6.5),
        legend.title=element_text(size=15))
Cycling_plot

#Endothelial cells
Breast_integrated$'Endothelial weight self-annotation'=newdata$Endothelial_weight
Endothelial_plot<-SpatialFeaturePlot(Breast_integrated,"Endothelial weight self-annotation",
                                 image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Endothelial_plot

#pDCs
Breast_integrated$'pDC weight self-annotation'=newdata$pDC_weight
pDC_plot<-SpatialFeaturePlot(Breast_integrated,"pDC weight self-annotation",
                                     image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=5.5),
        legend.title=element_text(size=15))
pDC_plot

#FOLPR2 positive macrophage weight
Breast_integrated$'FOLPR2 positive macrophage weight self-annotation'=newdata$FOLR2_positive_macrophage_weight
FOLPR2_pos<-SpatialFeaturePlot(Breast_integrated,'FOLPR2 positive macrophage weight self-annotation',
                                    image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
FOLPR2_pos

#FOLPR2 negative macrophage weight
Breast_integrated$'FOLPR2 negative macrophage weight self-annotation'=newdata$FOLR2_negative_macrophage_weight
FOLPR2_neg<-SpatialFeaturePlot(Breast_integrated,'FOLPR2 negative macrophage weight self-annotation',
                               image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
FOLPR2_neg

#Smooth muscle cell weight
Breast_integrated$'Smooth muscle cell weight self-annotation'=newdata$Smoothmuscle_weight
Smooth_plot<-SpatialFeaturePlot(Breast_integrated,'Smooth muscle cell weight self-annotation',
                               image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Smooth_plot

#CRYAB,SAA1 positive epithelial cells weight
Breast_integrated$'CRYAB,SAA1 positive epithelial cell weight self-annotation'=newdata$CRYAB_SAA1_positive_epithelial_weight
CRYAB_plot<-SpatialFeaturePlot(Breast_integrated,'CRYAB,SAA1 positive epithelial cell weight self-annotation',
                                image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
CRYAB_plot

#Cancer epithelial cells
Breast_integrated$'Cancer epithelial cell weight self-annotation'=newdata$Cancer_epithelial_cell_weight
Cancer_plot<-SpatialFeaturePlot(Breast_integrated,'Cancer epithelial cell weight self-annotation',
                               image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Cancer_plot

#Plasmablast weight
Breast_integrated$'Plasmablast weight self-annotation'=newdata$Plasmablasts_weight
Plasmablast_plot<-SpatialFeaturePlot(Breast_integrated,'Plasmablast weight self-annotation',
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Plasmablast_plot

#NK cell positive, T-cell positive weight
Breast_integrated$'NK cell positive, T-cell positive weight self-annotation'=newdata$NKcell_Tcell_weight
NK_plot<-SpatialFeaturePlot(Breast_integrated,'NK cell positive, T-cell positive weight self-annotation',
                                image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
NK_plot

#Myoepithelial weight
Breast_integrated$'Myoepithelial weight self-annotation'=newdata$Myoepithelial_weight
Myoepithelialplot<-SpatialFeaturePlot(Breast_integrated,"Myoepithelial weight self-annotation",
                                      image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Myoepithelialplot

#S100A8 positive epithelial cell weight
Breast_integrated$'S100A8 epithelial cells self-annotation'=newdata$S100A8_epithelial_cell_weight
S100A8_plot<-SpatialFeaturePlot(Breast_integrated,"S100A8 epithelial cells self-annotation",
                                      image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
S100A8_plot

#IGKV1-5 positive, T-cell rich weight
Breast_integrated$'IGKV1-5 positive, T-cell rich weight self-annotation'=newdata$IGKV15_Tcell_weight
IGKV1_plot<-SpatialFeaturePlot(Breast_integrated,"IGKV1-5 positive, T-cell rich weight self-annotation",
                                image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
IGKV1_plot

#SCGB1B2P positive epithelial cell weight
Breast_integrated$'SCGB1B2P positive epithelial cell weight self-annotation'=newdata$SCGB1B2P_epithelial_weight
SCG_plot<-SpatialFeaturePlot(Breast_integrated,"SCGB1B2P positive epithelial cell weight self-annotation",
                               image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
SCG_plot

####RCTD DECONVOLUTION OF BREAST TISSUE SAMPLE USING WU SINGLE-CELL ANNOTATIONS.
####Wu data was already loaded into the environment and pre-processed and libraries are loaded. Therefore these steps will be skipped. 
### Making single-cell references. It requires two compulsory inputs: 'counts' and 'cell_types'.
##Counts is a matrix representing the Digital Gene Expression
counts_1<-Wu.data
##Cell_types are the named factor of cell type for each cell.
#Retrieving meta data of Wu scRNA data (to retrieve author annotations)
meta_data_author<-read.csv(file.path(
  "~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/filtered_gene_bc_matrices/Wu_etal_2021_BRCA_scRNASeq","metadata.csv"))
#Obtaining cell type information 
cell_types_author<-meta_data_author$celltype_major
names(cell_types_author)<-meta_data_author$X
cell_type_author<-as.factor(cell_types_author)
#Constructing a single-cell reference
reference_author<-Reference(counts_1,cell_type_author)
###Since spatial transcriptomics and spatialRNA object is constructed, step is skipped. 
###Creating and running RCTD in full 'doubletmode' for 'major_cell_type' category within the Wu data. 
myRCTD_author<-create.RCTD(SpatialRNA,reference_author,max_cores=6)
myRCTD_author <- run.RCTD(myRCTD_author,doublet_mode = "full")
#Analyzing RCTD results
results_author<-myRCTD_author@results
#Normalizing cell type proportion to sum to 1
norm_weights_author=normalize_weights(results_author$weights)
#List of cell type names
cell_type_names_author<-myRCTD_author@cell_type_info$info[[2]]
spatialRNA_author<-myRCTD_author@spatialRNA
###RCTD results
##Cell contribution of spots mapped onto the original breast tissue biological image. 
#Making a dataframe with cell weights
newbarcode_author<-norm_weights_author@Dimnames[[1]]
nUMIs_author<-norm_weights_author@Dimnames[[1]]
newdata_author <- data.frame(Bcells_weight = norm_weights_author[,'B-cells'],CAFs_weight = norm_weights_author[, 'CAFs'],
                      Cancerepithelial_weight = norm_weights_author[, 'Cancer Epithelial'], 
                      Endothelial_weight = norm_weights_author[, 'Endothelial'],  
                      Myeloid_weight = norm_weights_author[, 'Myeloid'],  
                      NormalEpithelial_weight = norm_weights_author[, 'Normal Epithelial'], 
                      Plasmablast_weight = norm_weights_author[, 'Plasmablasts'], 
                      PVL_weight = norm_weights_author[, 'PVL'], 
                      Tcells_weight = norm_weights_author[, 'T-cells'])
#Add the new metadata to breast image
Breast_integrated_author <- AddMetaData(Breast, metadata = newdata_author)

#B-cells weight

Breast_integrated_author$'Bcell weight' = newdata_author$Bcells_weight
Bcellplot<-SpatialFeaturePlot(Breast_integrated_author,"Bcell weight",
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15),
  )
Bcellplot

#CAFs weight
Breast_integrated_author$'CAFs weight'=newdata_author$CAFs_weight
CAFplot<-SpatialFeaturePlot(Breast_integrated_author,"CAFs weight",
                            image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
CAFplot

#Cancer epithelial weight
Breast_integrated_author$'Cancer Epithelial weight'=newdata_author$Cancerepithelial_weight
Epithelialplot<-SpatialFeaturePlot(Breast_integrated_author,"Cancer Epithelial weight",
                                   image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Epithelialplot

#Endothelial weight
Breast_integrated_author$'Endothelial weight'=newdata_author$Endothelial_weight
Endothelialplot<-SpatialFeaturePlot(Breast_integrated_author,'Endothelial weight',
                                    image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Endothelialplot

#Myeloid weight
Breast_integrated_author$'Myeloid weight'=newdata_author$Myeloid_weight
Myeloidplot<-SpatialFeaturePlot(Breast_integrated_author,'Myeloid weight',
                                image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Myeloidplot

#Normal epithelial weight
Breast_integrated_author$'Normal Epithelial weight'=newdata_author$NormalEpithelial_weight
Normalepithelialplot<-SpatialFeaturePlot(Breast_integrated_author,'Normal Epithelial weight',
                                         image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Normalepithelialplot

#Plasmablasts weight
Breast_integrated_author$'Plasmablast weight'=newdata_author$Plasmablast_weight
Plasmablastplot<-SpatialFeaturePlot(Breast_integrated_author,'Plasmablast weight',image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Plasmablastplot

#PVL_weight
Breast_integrated_author$'PVL weight'=newdata_author$PVL_weight
PVLplot<-SpatialFeaturePlot(Breast_integrated_author,'PVL weight',
                            image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
PVLplot

#Tcells_weight
Breast_integrated_author$'Tcell weight'=newdata_author$Tcells_weight
Tcellplot<-SpatialFeaturePlot(Breast_integrated_author,'Tcell weight',
                              image.alpha = 0.25,alpha = c(1,20)) &
  theme(text = element_text(face = "bold"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=15))
Tcellplot





















































