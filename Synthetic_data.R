#######IMPORTING WU DATASET, DATA QUALITY CONTROL AND RANDOM DIVISION OF DATA FOR DECONVOLUTION INPUT AND SYNTHETIC DATA 

#####Loading necessary libraries.
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(synthspot)

####Loading single-cell data collected by Wu et al (2021),quality control of data and subsetting for synthetic spot production. 
Wu.data<-Read10X("~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/filtered_gene_bc_matrices/Wu_etal_2021_BRCA_scRNASeq",
                 gene.column = 1)
meta_data<-read.csv(file.path(
  "~/Desktop/Pharmacology_honours/DRT/Diss_milestone /Practice_tutorial/filtered_gene_bc_matrices/Wu_etal_2021_BRCA_scRNASeq","metadata.csv"))
Wu<-CreateSeuratObject(counts=Wu.data, project= "Wu3k",min.cells=3,min.features=200, meta.data = meta_data)
Wu[["percent.mito"]] <- PercentageFeatureSet(Wu,pattern = "^MT-")
###Filtering out poor data, indicative of cell stress or artifacts from cell-sorting process.
Wu<-subset(Wu,subset=nFeature_RNA>200 & nFeature_RNA<4000 & percent.mito<5)
###Sub-setting Wu data based on cancer sub-type.
Wu_HER2<-subset(Wu,subset=subtype=="HER2+")
Wu_ER<-subset(Wu,subset=subtype=="ER+")
Wu_TNBC<-subset(Wu,subset=subtype=="TNBC")

####Halving Wu collected cell signatures for cancer sub-types into 2 equal data sets. 1 for deconvolution tool input the other for synthspot for synthetic spot production

###Halving HER2 data into 2 equal subsets
synthspotCells <- sample(colnames(Wu_HER2),size = ncol(Wu_HER2)/2, replace = FALSE)
deconCells <- setdiff(colnames(Wu_HER2),synthspotCells)
#Create two Seurat object with the cells selected for synthspot and those kept for deconvolution
Wu_HER2_synthspot <- subset(Wu_HER2,subset = X %in% synthspotCells)
Wu_HER2_deconvolution <- subset(Wu_HER2, subset = X %in% deconCells)

###Halving ER data into 2 equal subsets
synthspotCells_ER <- sample(colnames(Wu_ER),size = ncol(Wu_ER)/2, replace = FALSE)
deconCells_ER <- setdiff(colnames(Wu_ER),synthspotCells_ER)
Wu_ER_synthspot <- subset(Wu_ER,subset = X %in% synthspotCells_ER)
Wu_ER_deconvolution <- subset(Wu_ER, subset = X %in% deconCells_ER)

###Halving TNBC data into 2 equal subsets
synthspotCells_TNBC <- sample(colnames(Wu_TNBC),size = ncol(Wu_TNBC)/2, replace = FALSE)
deconCells_TNBC <- setdiff(colnames(Wu_TNBC),synthspotCells_TNBC)
Wu_TNBC_synthspot <- subset(Wu_TNBC,subset = X %in% synthspotCells_TNBC)
Wu_TNBC_deconvolution <- subset(Wu_TNBC, subset = X %in% deconCells_TNBC)

####Processing both data subsets according to Seurat workflow (see single-cell annotation workflow).

###Processing Wu_HER data-sets

##Processing Wu_HER synthspot data
Wu_HER2_synthspot<-NormalizeData(Wu_HER2_synthspot)
Wu_HER2_synthspot<-FindVariableFeatures(Wu_HER2_synthspot,selection.method = "vst",nfeatures = 2000)
all.genes.Wu_HER2_synthspot<-rownames(Wu_HER2_synthspot)
Wu_HER2_synthspot<-ScaleData(Wu_HER2_synthspot,features=all.genes.Wu_HER2_synthspot)
Wu_HER2_synthspot<-RunPCA(Wu_HER2_synthspot,features=VariableFeatures(object=Wu_HER2_synthspot))
ElbowPlot(Wu_HER2_synthspot)
Wu_HER2_synthspot<-FindNeighbors(Wu_HER2_synthspot,dims=1:15)
Wu_HER2_synthspot<-FindClusters(Wu_HER2_synthspot,resolution = 0.10)
Wu_HER2_synthspot<-RunUMAP(Wu_HER2_synthspot,dims=1:15)
DimPlot(Wu_HER2_synthspot,reduction="umap")
Wu_HER2_synthspot[["RNA"]] <- as(Wu_HER2_synthspot[["RNA"]], Class = "Assay")

##Processing Wu_HER deconvolution data
Wu_HER2_deconvolution<-NormalizeData(Wu_HER2_deconvolution)
Wu_HER2_deconvolution<-FindVariableFeatures(Wu_HER2_deconvolution,selection.method = "vst",nfeatures = 2000)
all.genes.Wu_HER2_deconvolution<-rownames(Wu_HER2_deconvolution)
Wu_HER2_deconvolution<-ScaleData(Wu_HER2_deconvolution,features=all.genes.Wu_HER2_deconvolution)
Wu_HER2_deconvolution<-RunPCA(Wu_HER2_deconvolution,features=VariableFeatures(object=Wu_HER2_deconvolution))
ElbowPlot(Wu_HER2_deconvolution)
Wu_HER2_deconvolution<-FindNeighbors(Wu_HER2_deconvolution,dims=1:15)
WWu_HER2_deconvolution<-FindClusters(Wu_HER2_deconvolution,resolution = 0.10)
Wu_HER2_deconvolution<-RunUMAP(Wu_HER2_deconvolution,dims=1:15)
DimPlot(Wu_HER2_deconvolution,reduction="umap")
Wu_HER2_deconvolution[["RNA"]] <- as(Wu_HER2_deconvolution[["RNA"]], Class = "Assay")

###Processing Wu_ER data-sets. 

##Processing Wu_ER synthspot data
Wu_ER_synthspot<-NormalizeData(Wu_ER_synthspot)
Wu_ER_synthspot<-FindVariableFeatures(Wu_ER_synthspot,selection.method = "vst",nfeatures = 2000)
all.genes.Wu_ER_synthspot<-rownames(Wu_ER_synthspot)
Wu_ER_synthspot<-ScaleData(Wu_ER_synthspot,features=all.genes.Wu_ER_synthspot)
Wu_ER_synthspot<-RunPCA(Wu_ER_synthspot,features=VariableFeatures(object=Wu_ER_synthspot))
ElbowPlot(Wu_ER_synthspot)
Wu_ER_synthspot<-FindNeighbors(Wu_ER_synthspot,dims=1:15)
Wu_ER_synthspot<-FindClusters(Wu_ER_synthspot,resolution = 0.10)
Wu_ER_synthspot<-RunUMAP(Wu_ER_synthspot,dims=1:15)
DimPlot(Wu_ER_synthspot,reduction="umap")
Wu_ER_synthspot[["RNA"]] <- as(Wu_ER_synthspot[["RNA"]], Class = "Assay")

##Processing Wu_ER deconvolution data
Wu_ER_deconvolution<-NormalizeData(Wu_ER_deconvolution)
Wu_ER_deconvolution<-FindVariableFeatures(Wu_ER_deconvolution,selection.method = "vst",nfeatures = 2000)
all.genes.Wu_ER_deconvolution<-rownames(Wu_ER_deconvolution)
Wu_ER_deconvolution<-ScaleData(Wu_ER_deconvolution,features=all.genes.Wu_ER_deconvolution)
Wu_ER_deconvolution<-RunPCA(Wu_ER_deconvolution,features=VariableFeatures(object=Wu_ER_deconvolution))
ElbowPlot(Wu_ER_deconvolution)
Wu_ER_deconvolution<-FindNeighbors(Wu_ER_deconvolution,dims=1:15)
Wu_ER_deconvolution<-FindClusters(Wu_ER_deconvolution,resolution = 0.10)
Wu_ER_deconvolution<-RunUMAP(Wu_ER_deconvolution,dims=1:15)
DimPlot(Wu_ER_deconvolution,reduction="umap")
Wu_ER_deconvolution[["RNA"]] <- as(Wu_ER_deconvolution[["RNA"]], Class = "Assay")

###Processing Wu_TNBC data-sets

##Processing Wu_TNBC synthspot data
Wu_TNBC_synthspot<-NormalizeData(Wu_TNBC_synthspot)
Wu_TNBC_synthspot<-FindVariableFeatures(Wu_TNBC_synthspot,selection.method = "vst",nfeatures = 2000)
all.genes.Wu_TNBC_synthspot<-rownames(Wu_TNBC_synthspot)
Wu_TNBC_synthspot<-ScaleData(Wu_TNBC_synthspot,features=all.genes.Wu_TNBC_synthspot)
Wu_TNBC_synthspot<-RunPCA(Wu_TNBC_synthspot,features=VariableFeatures(object=Wu_TNBC_synthspot))
ElbowPlot(Wu_TNBC_synthspot)
Wu_TNBC_synthspot<-FindNeighbors(Wu_TNBC_synthspot,dims=1:15)
Wu_TNBC_synthspot<-FindClusters(Wu_TNBC_synthspot,resolution = 0.10)
Wu_TNBC_synthspot<-RunUMAP(Wu_TNBC_synthspot,dims=1:15)
DimPlot(Wu_TNBC_synthspot,reduction="umap")
Wu_TNBC_synthspot[["RNA"]] <- as(Wu_TNBC_synthspot[["RNA"]], Class = "Assay")

##Processing Wu_TNBC deconvolution data
Wu_TNBC_deconvolution<-NormalizeData(Wu_TNBC_deconvolution)
Wu_TNBC_deconvolution<-FindVariableFeatures(Wu_TNBC_deconvolution,selection.method = "vst",nfeatures = 2000)
all.genes.Wu_TNBC_deconvolution<-rownames(Wu_TNBC_deconvolution)
Wu_TNBC_deconvolution<-ScaleData(Wu_TNBC_deconvolution,features=all.genes.Wu_TNBC_deconvolution)
Wu_TNBC_deconvolution<-RunPCA(Wu_TNBC_deconvolution,features=VariableFeatures(object=Wu_TNBC_deconvolution))
ElbowPlot(Wu_TNBC_deconvolution)
Wu_TNBC_deconvolution<-FindNeighbors(Wu_TNBC_deconvolution,dims=1:15)
Wu_TNBC_deconvolution<-FindClusters(Wu_TNBC_deconvolution,resolution = 0.10)
Wu_TNBC_deconvolution<-RunUMAP(Wu_TNBC_deconvolution,dims=1:15)
DimPlot(Wu_TNBC_deconvolution,reduction="umap")
Wu_TNBC_deconvolution[["RNA"]] <- as(Wu_TNBC_deconvolution[["RNA"]], Class = "Assay")

####Generating synthetic data for each cancer subtype using synthspot data-sets. 
####Synthetic data for each cancer subtype is produced at 3 levels of resolution:major cell type, minor cell type, cell subset

###HER2+ synthetic data generation 

##HER2+ major cell type simulation 

#Producing a function to determine the dominant cell type in the synthetic HER2+ major cell type. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type. 
determine_dominant_celltype <- function(synthetic_visium_Wu_HER2_major) {
  spot_comp_HER_major<-synthetic_visium_Wu_HER2_major$spot_composition
  spot_comp_HER_major<-spot_comp_HER_major %>% 
    pivot_longer(cols = c(1:9),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts<-aggregate(Count~Cell_types,data=spot_comp_HER_major,FUN=sum)
  max_count_index <- which.max(total_counts$Count)
  dominant_celltype <- total_counts$Cell_types[max_count_index]
  print(dominant_celltype)}
#Generating synthetic data for HER2+ major cell type. The generation loops until the dominant cell type is 'Cancer.Epithelial'.
#The loop is constructed using the function above. 
synthetic_visium_Wu_HER2_major <- generate_synthetic_visium(seurat_obj = Wu_HER2_synthspot,
                                                            dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                            clust_var = "celltype_major",
                                                            n_regions = 5,
                                                            n_spots_min = 50, n_spots_max = 500,
                                                            n_cells_min = 3, n_cells_max = 30,
                                                            visium_mean = 30000, visium_sd = 8000)
dominant_celltype <- determine_dominant_celltype(synthetic_visium_Wu_HER2_major)
while (dominant_celltype != "Cancer.Epithelial") {
  synthetic_visium_Wu_HER2_major <- generate_synthetic_visium(seurat_obj = Wu_HER2_synthspot,
                                                              dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                              clust_var = "celltype_major",
                                                              n_regions = 5,
                                                              n_spots_min = 50, n_spots_max = 500,
                                                              n_cells_min = 3, n_cells_max = 30,
                                                              visium_mean = 30000, visium_sd = 8000)
  dominant_celltype <- determine_dominant_celltype(synthetic_visium_Wu_HER2_major)
}
##HER2+ minor cell type simulation 
#Synthspot has the tendency of potentially 'over-picking' low abundant cells for spot composition. 
#To circumvent this bias minor cell types containing less than 25 cells are removed. To do this 
#cell count is inspected using the command line below
celltype_minor_HER2_counts <- table(Wu_HER2_synthspot@meta.data$celltype_minor)
print(celltype_minor_HER2_counts)
#The following cell types are removed: 'B cells naive', 'Cancer Basal SC', 'Cancer LumA SC',
#'Cancer LumB SC', 'Cycling PVL','Cycling_Myeloid', 'Lymphatic LYVE1'. 
Wu_HER2_synthspot_minor<-subset(Wu_HER2_synthspot,subset = celltype_minor %in% c("B cells Memory","CAFs MSC iCAF-like","CAFs myCAF-like",
                                                                                 "Cancer Cycling","Cancer Her2 SC","Cycling T-cells",
                                                                                 "DCs","Endothelial ACKR1","Endothelial CXCL12",
                                                                                 "Endothelial RGS5","Luminal Progenitors","Macrophage",
                                                                                 "Mature Luminal","Monocyte","Myoepithelial","NK cells",
                                                                                 "NKT cells","Plasmablasts","PVL Differentiated","PVL Immature",
                                                                                 "T cells CD4+","T cells CD8+"))
#Producing a function to determine the dominant cell type in the synthetic HER2+ minor cell type. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.
determine_dominant_celltype_minor <- function(synthetic_visium_Wu_HER2_minor) {
  spot_comp_HER_minor<-synthetic_visium_Wu_HER2_minor$spot_composition
  spot_comp_HER_minor<-spot_comp_HER_minor %>% 
    pivot_longer(cols = c(1:22),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_minor<-aggregate(Count~Cell_types,data=spot_comp_HER_minor,FUN=sum)
  max_count_index_minor <- which.max(total_counts_minor$Count)
  dominant_celltype_minor <- total_counts_minor$Cell_types[max_count_index_minor]
  print(dominant_celltype_minor)
}
#Generating synthetic data for HER2+ minor cell type. The generation loops until the dominant cell type is 'Cancer.Her2.SC'.
#The loop is constructed using the function above. 
synthetic_visium_Wu_HER2_minor <- generate_synthetic_visium(seurat_obj = Wu_HER2_synthspot_minor,
                                                            dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                            clust_var = "celltype_minor",
                                                            n_regions = 5,
                                                            n_spots_min = 50, n_spots_max = 500,
                                                            n_cells_min = 3, n_cells_max = 30,
                                                            visium_mean = 30000, visium_sd = 8000)
dominant_celltype_minor <- determine_dominant_celltype_minor(synthetic_visium_Wu_HER2_minor)
while (dominant_celltype_minor != "Cancer.Her2.SC") {
  synthetic_visium_Wu_HER2_minor <- generate_synthetic_visium(seurat_obj = Wu_HER2_synthspot_minor,
                                                              dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                              clust_var = "celltype_minor",
                                                              n_regions = 5,
                                                              n_spots_min = 50, n_spots_max = 500,
                                                              n_cells_min = 3, n_cells_max = 30,
                                                              visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_minor <- determine_dominant_celltype_minor(synthetic_visium_Wu_HER2_minor)
}

##HER2+ subset cell type simulation
#Inspecting counts within each cell subtype class found in the HER2+ subset data.
celltype_subset_HER2_counts <- table(Wu_HER2_synthspot@meta.data$celltype_subset)
print(celltype_subset_HER2_counts)
#The following cell subsets contain less than 25 cells: 'B cells Naive', 'Cancer Basal SC','Cancer LumA SC',
#'Cancer LumB SC','Cycling PVL','Cycling_Myeloid','Endothelial Lymphatic LYVE1','Myeloid_c0_DC_LAMP3',
#'Myeloid_c11_cDC2_CD1C','Myeloid_c3_cDC1_CLEC9A','Myeloid_c9_Macrophage_2_CXCL10','PVL_Immature s2',
#'T_cells_c5_CD8+_GZMK'. 
#Cells aforementioned are removed via the function below. 
Wu_HER2_synthspot_subset<-subset(Wu_HER2_synthspot,subset=celltype_subset %in% c("B cells Memory",
                        "CAFs MSC iCAF-like s1","CAFs MSC iCAF-like s2","CAFs myCAF like s4","CAFs myCAF like s5",
                        "CAFs Transitioning s3","Cancer Cycling","Cancer Her2 SC","Endothelial ACKR1",
                        "Endothelial CXCL12","Endothelial RGS5","Luminal Progenitors","Mature Luminal",
                        "Myeloid_c1_LAM1_FABP5","Myeloid_c10_Macrophage_1_EGR1","Myeloid_c12_Monocyte_1_IL1B",
                        "Myeloid_c2_LAM2_APOE","Myeloid_c4_DCs_pDC_IRF7","Myeloid_c8_Monocyte_2_S100A9",
                        "Myoepithelial","Plasmablasts","PVL Differentiated s3",
                        "PVL Immature s1","T_cells_c0_CD4+_CCR7","T_cells_c1_CD4+_IL7R",
                        "T_cells_c10_NKT_cells_FCGR3A","T_cells_c11_MKI67","T_cells_c2_CD4+_T-regs_FOXP3",
                        "T_cells_c3_CD4+_Tfh_CXCL13","T_cells_c4_CD8+_ZFP36","T_cells_c6_IFIT1",
                        "T_cells_c7_CD8+_IFNG","T_cells_c8_CD8+_LAG3","T_cells_c9_NK_cells_AREG"))

#Producing a function to determine the dominant cell type in the synthetic HER2+ cell subset data. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.
determine_dominant_celltype_subset <- function(synthetic_visium_Wu_HER2_subset) {
  spot_comp_HER_subset<-synthetic_visium_Wu_HER2_subset$spot_composition
  spot_comp_HER_subset<-spot_comp_HER_subset %>% 
    pivot_longer(cols = c(1:34),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_subset<-aggregate(Count~Cell_types,data=spot_comp_HER_subset,FUN=sum)
  max_count_index_subset <- which.max(total_counts_subset$Count)
  dominant_celltype_subset <- total_counts_subset$Cell_types[max_count_index_subset]
  print(dominant_celltype_subset)
}
#Generating synthetic data for HER2+ cell subset type. The generation loops until the dominant cell type is 'Cancer.Her2.SC'.
#The loop is constructed using the function above.  
synthetic_visium_Wu_HER2_subset <- generate_synthetic_visium(seurat_obj = Wu_HER2_synthspot_subset,
                                                            dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                            clust_var = "celltype_subset",
                                                            n_regions = 5,
                                                            n_spots_min = 50, n_spots_max = 500,
                                                            n_cells_min = 3, n_cells_max = 30,
                                                            visium_mean = 30000, visium_sd = 8000)
dominant_celltype_subset <- determine_dominant_celltype_subset(synthetic_visium_Wu_HER2_subset)
while (dominant_celltype_subset != "Cancer.Her2.SC") {
  synthetic_visium_Wu_HER2_subset <- generate_synthetic_visium(seurat_obj = Wu_HER2_synthspot_subset,
                                                              dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                              clust_var = "celltype_subset",
                                                              n_regions = 5,
                                                              n_spots_min = 50, n_spots_max = 500,
                                                              n_cells_min = 3, n_cells_max = 30,
                                                              visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_subset <- determine_dominant_celltype_subset(synthetic_visium_Wu_HER2_subset)
}

###ER+ synthetic data generation.

##ER2+ major cell type simulation
#Producing a function to determine the dominant cell type in the synthetic ER+ major cell type subset. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.
determine_dominant_celltype_ERmajor <- function(synthetic_visium_Wu_ER_major) {
  spot_comp_ER_major<-synthetic_visium_Wu_ER_major$spot_composition
  spot_comp_ER_major<-spot_comp_ER_major %>% 
    pivot_longer(cols = c(1:9),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_ER_major<-aggregate(Count~Cell_types,data=spot_comp_ER_major,FUN=sum)
  max_count_index_ER_major <- which.max(total_counts_ER_major$Count)
  dominant_celltype_ER_major <- total_counts_ER_major$Cell_types[max_count_index_ER_major]
  print(dominant_celltype_ER_major)
}
#Generating synthetic data for ER+ major cell type. The generation loops until the dominant cell type is 'Cancer.Epithelial'.
#The loop is constructed using the function above.
synthetic_visium_Wu_ER_major <- generate_synthetic_visium(seurat_obj = Wu_ER_synthspot,
                                                             dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                             clust_var = "celltype_major",
                                                             n_regions = 5,
                                                             n_spots_min = 50, n_spots_max = 500,
                                                             n_cells_min = 3, n_cells_max = 30,
                                                             visium_mean = 30000, visium_sd = 8000)
dominant_celltype_ER_major <- determine_dominant_celltype_ERmajor(synthetic_visium_Wu_ER_major)
while (dominant_celltype_ER_major != "Cancer.Epithelial") {
  synthetic_visium_Wu_ER_major <- generate_synthetic_visium(seurat_obj = Wu_ER_synthspot,
                                                               dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                               clust_var = "celltype_major",
                                                               n_regions = 5,
                                                               n_spots_min = 50, n_spots_max = 500,
                                                               n_cells_min = 3, n_cells_max = 30,
                                                               visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_ER_major <- determine_dominant_celltype_ERmajor(synthetic_visium_Wu_ER_major)
}

##ER+ minor cell type simulation
#Inspecting counts within each cell subtype class found in the ER+ major data.
celltype_minor_ER_counts <- table(Wu_ER_synthspot@meta.data$celltype_minor)
print(celltype_minor_ER_counts)
#The following cell subsets contain less than 25 cells: 'Cancer Basal SC','Cancer Her2 SC','Cycling PVL',
#'Cycling_Myeloid','Endothelial Lymphatic LYVE1'.
#Cells aforementioned are removed via the function below. 
Wu_ER_synthspot_minor<-subset(Wu_ER_synthspot,subset=celltype_minor %in% c("B cells Memory","CAFs MSC iCAF-like","CAFs myCAF-like",
                                                                           "Cancer Cycling","Cancer LumA SC","Cancer LumB SC",
                                                                           "Cycling T-cells","DCs","Endothelial ACKR1","Endothelial CXCL12",
                                                                           "Endothelial RGS5","Luminal Progenitors","Macrophage",
                                                                           "Mature Luminal","Monocyte","Myoepithelial","NK cells",
                                                                           "NKT cells","Plasmablasts","PVL Differentiated","PVL immature",
                                                                           "T cells CD4+","T cells CD8+"))
#Producing a function to determine the dominant cell type in the synthetic ER+ minor cell type subset. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.
determine_dominant_celltype_ERminor <- function(synthetic_visium_Wu_ER_minor) {
  spot_comp_ER_minor<-synthetic_visium_Wu_ER_minor$spot_composition
  spot_comp_ER_minor<-spot_comp_ER_minor %>% 
    pivot_longer(cols = c(1:22),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_ER_minor<-aggregate(Count~Cell_types,data=spot_comp_ER_minor,FUN=sum)
  max_count_index_ER_minor<- which.max(total_counts_ER_minor$Count)
  dominant_celltype_ER_minor <- total_counts_ER_minor$Cell_types[max_count_index_ER_minor]
  print(dominant_celltype_ER_minor)
}

#Generating synthetic data for ER+ minor cell type. The generation loops until the dominant cell type is 'Cancer.LumA.SC'.
#The loop is constructed using the function above.
synthetic_visium_Wu_ER_minor <- generate_synthetic_visium(seurat_obj = Wu_ER_synthspot_minor,
                                                          dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                          clust_var = "celltype_minor",
                                                          n_regions = 5,
                                                          n_spots_min = 50, n_spots_max = 500,
                                                          n_cells_min = 3, n_cells_max = 30,
                                                          visium_mean = 30000, visium_sd = 8000)
dominant_celltype_ER_minor <- determine_dominant_celltype_ERminor(synthetic_visium_Wu_ER_minor)
while (dominant_celltype_ER_minor != "Cancer.LumA.SC") {
  synthetic_visium_Wu_ER_minor <- generate_synthetic_visium(seurat_obj = Wu_ER_synthspot_minor,
                                                              dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                              clust_var = "celltype_minor",
                                                              n_regions = 5,
                                                              n_spots_min = 50, n_spots_max = 500,
                                                              n_cells_min = 3, n_cells_max = 30,
                                                              visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_ER_minor <- determine_dominant_celltype_ERminor(synthetic_visium_Wu_ER_minor)
}

##ER+ subset cell type simulation 
#Inspecting counts within each cell subtype class found in the ER data
celltype_subset_ER_counts <- table(Wu_ER_synthspot@meta.data$celltype_subset)
print(celltype_subset_ER_counts)
#Removing cell sub types that are under 25 due to previously mentioned reason. Cells removed from dataset: 
#'CAFs MSC iCAF-like s2', 'CAFs myCAF like s4', 'Cancer Basal SC', 'Cancer HER2 SC', 'Cycling PVL', 
#'Cycling_Myeloid', 'Endothelial Lymphatic LYVE1', 'Myeloid_c0_DC-LAMP3', 'Myeloid_c3_cDC1_CLEC9A',
#'Myleoid_c7_Monocyte_3_FCGR3A', 'Myeloid_c9_Macrophage_2_CXCL10'. 
Wu_ER_synthspot_subset<-subset(Wu_ER_synthspot,subset=celltype_subset %in% c("B cells Memory","CAFs MSC iCAF-like s1", "CAFs myCAF like s5",
                                                                             "CAFs transitioning s3","Cancer Cycling","Cancer LumA SC","Cancer LumB SC","Endothelial ACKR1",
                                                                             "Endothelial CXCL12","Endothelial RGS5","Luminal Progenitors","Mature Luminal",
                                                                             "Myeloid_c1_LAM1_FABP5","Myeloid_c10_Macrophage_1_ERG1","Myeloid_c11_cDC2_CD1C",
                                                                             "Myeloid_c12_Monocyte_1_IL1B","Myeloid_c2_LAM2_APOE","Myeloid_c4_DCs_pDC_IRF7",
                                                                             "Myeloid_c8_Monocyte_2_S100A9","Myoepithelial","Plasmablasts","PVL Differentiated s3",
                                                                             "PVL Immature s1","PVL_Immature s2","T_cells_c0_CD4+_CCR7","T_cells_c1_CD4+_IL7R",
                                                                             "T_cells_c10_NKT_cells_FCGR3A","T_cells_c11_MKI67","T_cells_c2_CD4+_T-regs_FOXP3",
                                                                             "T_cells_c3_CD4+_Tfh_CXCL13","T_cells_c4_CD8+_ZFP36","T_cells_c6_IFIT1",
                                                                             "T_cells_c7_CD8+_IFNG","T_cells_c8_CD8+_LAG3","T_cells_c9_NK_cells_AREG"))
#Producing a function to determine the dominant cell type in the synthetic ER+ cell subset data. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.
determine_dominant_celltype_subset_ER <- function(synthetic_visium_Wu_ER_subset) {
  spot_comp_ER_subset<-synthetic_visium_Wu_ER_subset$spot_composition
  spot_comp_ER_subset<-spot_comp_ER_subset %>% 
    pivot_longer(cols = c(1:33),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_subset_ER<-aggregate(Count~Cell_types,data=spot_comp_ER_subset,FUN=sum)
  max_count_index_subset_ER <- which.max(total_counts_subset_ER$Count)
  dominant_celltype_subset_ER <- total_counts_subset_ER$Cell_types[max_count_index_subset_ER]
  print(dominant_celltype_subset_ER)
}
#Generating synthetic data for ER+ cell subset data. The generation loops until the dominant cell type is 'Cancer.LumA.SC'.
#The loop is constructed using the function above.
synthetic_visium_Wu_ER_subset <- generate_synthetic_visium(seurat_obj = Wu_ER_synthspot_subset,
                                                             dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                             clust_var = "celltype_subset",
                                                             n_regions = 5,
                                                             n_spots_min = 50, n_spots_max = 500,
                                                             n_cells_min = 3, n_cells_max = 30,
                                                             visium_mean = 30000, visium_sd = 8000)
dominant_celltype_subset_ER <- determine_dominant_celltype_subset_ER(synthetic_visium_Wu_ER_subset)
while (dominant_celltype_subset_ER != "Cancer.LumA.SC") {
  synthetic_visium_Wu_ER_subset <- generate_synthetic_visium(seurat_obj = Wu_ER_synthspot_subset,
                                                               dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                               clust_var = "celltype_subset",
                                                               n_regions = 5,
                                                               n_spots_min = 50, n_spots_max = 500,
                                                               n_cells_min = 3, n_cells_max = 30,
                                                               visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_subset_ER <- determine_dominant_celltype_subset_ER(synthetic_visium_Wu_ER_subset)
}

###TNBC synthetic data generation.
##TNBC major cell type simulation
#Producing a function to determine the dominant cell type in the synthetic TNBC major cell type. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.
determine_dominant_celltype_TNBCmajor <- function(synthetic_visium_Wu_TNBC_major) {
  spot_comp_TNBC_major<-synthetic_visium_Wu_TNBC_major$spot_composition
  spot_comp_TNBC_major<-spot_comp_TNBC_major %>% 
    pivot_longer(cols = c(1:9),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_TNBC_major<-aggregate(Count~Cell_types,data=spot_comp_TNBC_major,FUN=sum)
  max_count_index_TNBC_major <- which.max(total_counts_TNBC_major$Count)
  dominant_celltype_TNBC_major <- total_counts_TNBC_major$Cell_types[max_count_index_TNBC_major]
  print(dominant_celltype_TNBC_major)
}
#Generating synthetic data for TNBC major cell type. The generation loops until the dominant cell type is 'Cancer.Epithelial'.
#The loop is constructed using the function above.
synthetic_visium_Wu_TNBC_major <- generate_synthetic_visium(seurat_obj = Wu_TNBC_synthspot,
                                                          dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                          clust_var = "celltype_major",
                                                          n_regions = 5,
                                                          n_spots_min = 50, n_spots_max = 500,
                                                          n_cells_min = 3, n_cells_max = 30,
                                                          visium_mean = 30000, visium_sd = 8000)
dominant_celltype_TNBC_major <- determine_dominant_celltype_TNBCmajor(synthetic_visium_Wu_TNBC_major)
while (dominant_celltype_TNBC_major != "Cancer.Epithelial") {
  synthetic_visium_Wu_TNBC_major <- generate_synthetic_visium(seurat_obj = Wu_TNBC_synthspot,
                                                            dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                            clust_var = "celltype_major",
                                                            n_regions = 5,
                                                            n_spots_min = 50, n_spots_max = 500,
                                                            n_cells_min = 3, n_cells_max = 30,
                                                            visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_TNBC_major <- determine_dominant_celltype_TNBCmajor(synthetic_visium_Wu_TNBC_major)
}

##TNBC minor cell type simulation. 
#Inspecting counts within each cell subtype class found in the TNBC minor data.
celltype_minor_TNBC_counts <- table(Wu_TNBC_synthspot@meta.data$celltype_minor)
print(celltype_minor_TNBC_counts)
#The following cells minor types have been found to possess less than 25 cell counts: 'Cancer LumA SC','Cycling PVL',
#'Lymphatic LYVE1'. These are removed from the TNBC minor data set
Wu_TNBC_synthspot_minor<-subset(Wu_TNBC_synthspot,subset = celltype_minor %in% c("B cells Memory","B cells Naive","CAFs MSC iCAF-like",
                                                                                 "CAFs myCAF-like","Cancer Basal SC","Cancer Cycling",
                                                                                 "Cancer Her2 SC","Cancer LumB SC","Cycling PVL","Cycling T-cells",
                                                                                 "Cycling_Myeloid","DCs","Endothelial ACKR1","Endothelial CXCL12",
                                                                                 "Endothelial RGS5","Luminal Progenitors","Macrophage","Mature Luminal",
                                                                                 "Monocyte","Myoepithelial","NK cells","NKT cells","Plasmablasts",
                                                                                 "PVL Differentiated","PVL Immature","T cells CD4+","T cells CD8+"))
#Producing a function to determine the dominant cell type in the synthetic TNBC minor cell type. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.                                        
determine_dominant_celltype_TNBCminor <- function(synthetic_visium_Wu_TNBCminor) {
  spot_comp_TNBCminor<-synthetic_visium_Wu_TNBCminor$spot_composition
  spot_comp_TNBCminor<-spot_comp_TNBCminor %>% 
    pivot_longer(cols = c(1:27),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_TNBCminor<-aggregate(Count~Cell_types,data=spot_comp_TNBCminor,FUN=sum)
  max_count_index_TNBCminor <- which.max(total_counts_TNBCminor$Count)
  dominant_celltype_TNBCminor <- total_counts_TNBCminor$Cell_types[max_count_index_TNBCminor]
  print(dominant_celltype_TNBCminor)
}
#Generating synthetic data for TNBC minor cell type. The generation loops until the dominant cell type is 'Cancer.Basal.SC'.
#The loop is constructed using the function above.
synthetic_visium_Wu_TNBCminor <- generate_synthetic_visium(seurat_obj = Wu_TNBC_synthspot_minor,
                                                            dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                            clust_var = "celltype_minor",
                                                            n_regions = 5,
                                                            n_spots_min = 50, n_spots_max = 500,
                                                            n_cells_min = 3, n_cells_max = 30,
                                                            visium_mean = 30000, visium_sd = 8000)
dominant_celltype_TNBCminor <- determine_dominant_celltype_TNBCminor(synthetic_visium_Wu_TNBCminor)
while (dominant_celltype_TNBCminor != "Cancer.Basal.SC") {
  synthetic_visium_Wu_TNBCminor <- generate_synthetic_visium(seurat_obj = Wu_TNBC_synthspot_minor,
                                                              dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                              clust_var = "celltype_minor",
                                                              n_regions = 5,
                                                              n_spots_min = 50, n_spots_max = 500,
                                                              n_cells_min = 3, n_cells_max = 30,
                                                              visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_TNBCminor <- determine_dominant_celltype_TNBCminor(synthetic_visium_Wu_TNBCminor)
}

##TNBC subset cell type simulation. 
#Inspecting counts within each cell subtype class found in the TNBC subset data.
celltype_subset_TNBC_counts <- table(Wu_TNBC_synthspot@meta.data$celltype_subset)
print(celltype_subset_TNBC_counts)
#The following cell sub types were found to have a cell count under 25: 'Cancer LumA SC', 'Cancer LumB SC','Cycling PVL',
#'Endothelial Lymphatic LYVE1', 'Myeloid_c0_DC_LAMP3','Myeloid_c5_Macrophage_3_SIGLEC1','Myeloid_c7_Monocyte_3_FCGR3A'. 
#These are removed from the subset dataset. 
Wu_TNBC_synthspot_subset<-subset(Wu_TNBC_synthspot,subset = celltype_subset %in% c("B cells Memory","B cells Naive","CAFs MSC iCAF-like s1",
                                                                                   "CAFs MSC iCAF-like s2","CAFs myCAF like s4","CAFs myCAF like s5",
                                                                                   "CAFs Transitioning s3","Cancer Basal SC","Cancer Cycling","Cancer Her2 SC",
                                                                                   "Cycling_Myeloid","Endothelial ACKR1","Endothelial CXCL12","Endothelial RGS5",
                                                                                   "Luminal Progenitors","Mature Luminal","Myeloid_c1_LAM1_FABP5","Myeloid_c10_Macrophage_1_EGR1",
                                                                                   "Myeloid_c11_cDC2_CD1C","Myeloid_c12_Monocyte_1_IL1B","Myeloid_c2_LAM2_APOE",
                                                                                   "Myeloid_3_cDC1_CLEC9A","Myeloid_c4_DCs_pDC_IRF7","Myeloid_c8_Monocyte_2_S100A9",
                                                                                   "Myeloid_c9_Macrophage_2_CXCL10","Plasmablasts","Myoepithelial",
                                                                                   "PVL Differentiated s3","PVL Immature S1","PVL_Immature s2","T_cells_c0_CD4+_CCR7",
                                                                                   "T_cells_c1_CD4+_IL7R","T_cells_c10_NKT_cells_FCGR3A","T_cells_c11_MKI67","T_cells_c2_CD4+_T-regs_FOXP3",
                                                                                   "T_cells_c3_CD4+_Tfh_CXCL13","T_cells_c4_CD8+_ZFP36","T_cells_c5_CD8+_GZMK","T_cells_c6_IFIT1",
                                                                                   "T_cells_c7_CD8+_IFNG","T_cells_c8_CD8+_LAG3","T_cells_c9_NK_cells_AREG"))
#Producing a function to determine the dominant cell type in the synthetic TNBC subset cell type. Function first collates cell_type counts.
#Secondly it finds the index of the cell type with the largest count and determines it to be the dominant cell type.    
determine_dominant_celltype_TNBCsubset <- function(synthetic_visium_Wu_TNBCsubset) {
  #Collating cell_type count
  spot_comp_TNBCsubset<-synthetic_visium_Wu_TNBCsubset$spot_composition
  spot_comp_TNBCsubset<-spot_comp_TNBCsubset %>% 
    pivot_longer(cols = c(1:40),
                 names_to = "Cell_types",
                 values_to = "Count")
  total_counts_TNBCsubset<-aggregate(Count~Cell_types,data=spot_comp_TNBCsubset,FUN=sum)
  # Find the index of the cell type with the maximum count
  max_count_index_TNBCsubset <- which.max(total_counts_TNBCsubset$Count)
  # Extract the cell type with the maximum count
  dominant_celltype_TNBCsubset <- total_counts_TNBCsubset$Cell_types[max_count_index_TNBCsubset]
  print(dominant_celltype_TNBCsubset)
}
#Generating synthetic data for TNBC subset cell type. The generation loops until the dominant cell type is 'Cancer.Basal.SC'.
#The loop is constructed using the function above.
synthetic_visium_Wu_TNBCsubset <- generate_synthetic_visium(seurat_obj = Wu_TNBC_synthspot_subset,
                                                           dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                           clust_var = "celltype_subset",
                                                           n_regions = 5,
                                                           n_spots_min = 50, n_spots_max = 500,
                                                           n_cells_min = 3, n_cells_max = 30,
                                                           visium_mean = 30000, visium_sd = 8000)
dominant_celltype_TNBCsubset <- determine_dominant_celltype_TNBCsubset(synthetic_visium_Wu_TNBCsubset)
while (dominant_celltype_TNBCsubset != "Cancer.Basal.SC") {
  synthetic_visium_Wu_TNBCsubset <- generate_synthetic_visium(seurat_obj = Wu_TNBC_synthspot_subset,
                                                             dataset_type = "artificial_partially_dominant_celltype_diverse",
                                                             clust_var = "celltype_subset",
                                                             n_regions = 5,
                                                             n_spots_min = 50, n_spots_max = 500,
                                                             n_cells_min = 3, n_cells_max = 30,
                                                             visium_mean = 30000, visium_sd = 8000)
  dominant_celltype_TNBCsubset <- determine_dominant_celltype_TNBCsubset(synthetic_visium_Wu_TNBCsubset)
}

####Plotting cell composition of all synthetic cancer subsets at all different levels of resolution (major,minor and subset)        
###HER2 cancer
##HER2 major
Plot_HER2_major<-synthetic_visium_Wu_HER2_major$spot_composition
Plot_HER2_major$name<-as.factor(Plot_HER2_major$name)
Plot_HER2_major<-Plot_HER2_major %>% 
  pivot_longer(cols=c("B.cells","CAFs","Cancer.Epithelial","Endothelial","Myeloid",
                      "Normal.Epithelial","Plasmablasts","PVL","T.cells"),
               names_to = "Cell_type_major",
               values_to = "Count")
Plot_HER2_major<-Plot_HER2_major %>% 
  ggplot(aes(x=Cell_type_major,y=Count)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "HER2 major")
Plot_HER2_major
##HER2 minor
unique(Wu_HER2_synthspot_minor@meta.data$celltype_minor)
Plot_HER2_minor<-synthetic_visium_Wu_HER2_minor$spot_composition
Plot_HER2_minor<-Plot_HER2_minor%>% 
  pivot_longer(cols = c(1:22),
               names_to = "Cell_types_minor",
               values_to = "Count")
Plot_HER2_minor<-Plot_HER2_minor %>% 
  ggplot(aes(x=Cell_types_minor,y=Count)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "HER2 minor")
Plot_HER2_minor
##HER2 subset
Plot_HER2_subset<-synthetic_visium_Wu_HER2_subset$spot_composition
Plot_HER2_subset<-Plot_HER2_subset%>% 
  pivot_longer(cols = c(1:32),
               names_to = "Cell_types_subtype",
               values_to = "Count")
Plot_HER2_subset<-Plot_HER2_subset %>% 
  ggplot(aes(x=Cell_types_subtype,y=Count)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "HER2 subset")
Plot_HER2_subset
###ER cancer
##ER major
Plot_ER_major<-synthetic_visium_Wu_ER_major$spot_composition
Plot_ER_major$name<-as.factor(Plot_ER_major$name)
Plot_ER_major<-Plot_ER_major %>% 
  pivot_longer(cols=c("B.cells","CAFs","Cancer.Epithelial","Endothelial","Myeloid",
                      "Normal.Epithelial","Plasmablasts","PVL","T.cells"),
               names_to = "Cell_type_major",
               values_to = "Count")
Plot_ER_major<-Plot_ER_major %>% 
  ggplot(aes(x=Cell_type_major,y=Count)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "ER major")
Plot_ER_major
##ER minor
Plot_ER_minor<-synthetic_visium_Wu_ER_minor$spot_composition
Plot_ER_minor<-Plot_ER_minor%>% 
  pivot_longer(cols = c(1:22),
               names_to = "Cell_types_minor",
               values_to = "Count")
Plot_ER_minor<-Plot_ER_minor %>% 
  ggplot(aes(x=Cell_types_minor,y=Count)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "ER minor")
Plot_ER_minor
##ER subset
Plot_ER_subset<-synthetic_visium_Wu_ER_subset$spot_composition
Plot_ER_subset<-Plot_ER_subset%>% 
  pivot_longer(cols = c(1:33),
               names_to = "Cell_types_subtype",
               values_to = "Count")
Plot_ER_subset<-Plot_ER_subset %>% 
  ggplot(aes(x=Cell_types_subtype,y=Count)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "ER subset")
Plot_ER_subset
###TNBC
##TNBC major
Plot_TNBC_major<-synthetic_visium_Wu_TNBC_major$spot_composition
Plot_TNBC_major$name<-as.factor(Plot_TNBC_major$name)
Plot_TNBC_major<-Plot_TNBC_major %>% 
  pivot_longer(cols=c("B.cells","CAFs","Cancer.Epithelial","Endothelial","Myeloid",
                      "Normal.Epithelial","Plasmablasts","PVL","T.cells"),
               names_to = "Cell_type_major",
               values_to = "Count")
Plot_TNBC_major<-Plot_TNBC_major %>% 
  ggplot(aes(x=Cell_type_major,y=Count)) +
  geom_boxplot() +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "TNBC major")
Plot_TNBC_major
##TNBC minor
Plot_TNBC_minor<-synthetic_visium_Wu_TNBCminor$spot_composition
Plot_TNBC_minor<-Plot_TNBC_minor%>% 
  pivot_longer(cols = c(1:27),
               names_to = "Cell_types_minor",
               values_to = "Count")
Plot_TNBC_minor<-Plot_TNBC_minor %>% 
  ggplot(aes(x=Cell_types_minor,y=Count)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "TNBC minor")
Plot_TNBC_minor
##TNBC subset
Plot_TNBC_subset<-synthetic_visium_Wu_TNBCsubset$spot_composition
Plot_TNBC_subset<-Plot_TNBC_subset%>% 
  pivot_longer(cols = c(1:40),
               names_to = "Cell_types_subtype",
               values_to = "Count")
Plot_TNBC_subset<-Plot_TNBC_subset %>% 
  ggplot(aes(x=Cell_types_subtype,y=Count)) +
  geom_boxplot() + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "TNBC subset")
Plot_TNBC_subset

###Saving synthetic data and deconvolution

##Saving deconvolution inputs
save(Wu_HER2_deconvolution,file="HER2_Deconvolution")
save(Wu_ER_deconvolution,file="ER_Deconvolution")
save(Wu_TNBC_deconvolution,file="TNBC_Deconvolution")
##Saving synthetic data sets produced
save(synthetic_visium_Wu_ER_major,file="synthetic_ER_major")
save(synthetic_visium_Wu_ER_minor,file="synthetic_ER_minor")
save(synthetic_visium_Wu_ER_subset,file="synthetic_ER_subset")
save(synthetic_visium_Wu_HER2_major,file="synthetic_HER2_major")
save(synthetic_visium_Wu_HER2_minor,file="synthetic_HER2_minor")
save(synthetic_visium_Wu_HER2_subset,file="synthetic_HER2_subset")
save(synthetic_visium_Wu_TNBC_major,file="synthetic_TNBC_major")
save(synthetic_visium_Wu_TNBCminor,file="synthetic_TNBC_minor")
save(synthetic_visium_Wu_TNBCsubset,file="synthetic_TNBC_subset")


