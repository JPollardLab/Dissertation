##WORKFLOW FOR OBTAINING RESULTS FOR: 'MAINTENANCE OF STABILITY FOR MORE COMPLEX PHENOTYPIC SUBSETS IN TNBC'
#Merging datasets obtained from running RCTD on the major and subset resolution level.
merged_df<-merge(RCTD_major,RCTD_subset,by="spotID",suffixes = c("_major","_subset"))

##B_cell resolution calculation
merged_df$B_cell_res<-abs(merged_df$rctd_B.cells - (merged_df$rctd_B.cells.Memory + merged_df$rctd_B.cells.Naive))
#Median B_cell resolution calculation
Bcell_stability<-median(merged_df$B_cell_res)

##CAFs resolution calculation
merged_df$CAF_res<-abs(merged_df$rctd_CAFs-(merged_df$rctd_CAFs.MSC.iCAF.like.s1+merged_df$rctd_CAFs.MSC.iCAF.like.s2
                                            +merged_df$rctd_CAFs.myCAF.like.s4+merged_df$rctd_CAFs.myCAF.like.s5+
                                              merged_df$rctd_CAFs.Transitioning.s3))
#Median CAF resolution calculation
CAF_stability<-median(merged_df$CAF_res)

##Cancer epithelial resolution calculation
merged_df$Cancer_res<-abs(merged_df$rctd_Cancer.Epithelial-(merged_df$rctd_Cancer.Basal.SC+
                                                                          merged_df$rctd_Cancer.Cycling+
                                                                          merged_df$rctd_Cancer.Her2.SC+
                                                                          merged_df$rctd_Cancer.LumA.SC+
                                                                          merged_df$rctd_Cancer.LumB.SC))
#Median cancer epithelial calculation
Cancer_Stability<-median(merged_df$Cancer_res)

##Endothelial resolution calculation
merged_df$Endothelial_res<-abs(merged_df$rctd_Endothelial-(merged_df$rctd_Endothelial.ACKR1+
                                                             merged_df$rctd_Endothelial.CXCL12+
                                                             merged_df$rctd_Endothelial.Lymphatic.LYVE1+
                                                             merged_df$rctd_Endothelial.RGS5))
#Median endothelial calculation
Endothelial_stability<-median(merged_df$Endothelial_res)

##Myeloid resolution calculation
merged_df$Myeloid_res<-abs(merged_df$rctd_Myeloid-(merged_df$rctd_Myeloid_c0_DC_LAMP3+
                                                     merged_df$rctd_Myeloid_c10_Macrophage_1_EGR1+
                                                     merged_df$rctd_Myeloid_c11_cDC2_CD1C+
                                                     merged_df$rctd_Myeloid_c12_Monocyte_1_IL1B+
                                                     merged_df$rctd_Myeloid_c1_LAM1_FABP5+
                                                     merged_df$rctd_Myeloid_c2_LAM2_APOE+
                                                     merged_df$rctd_Myeloid_c3_cDC1_CLEC9A+
                                                     merged_df$rctd_Myeloid_c4_DCs_pDC_IRF7+
                                                     merged_df$rctd_Myeloid_c5_Macrophage_3_SIGLEC1+
                                                     merged_df$rctd_Myeloid_c7_Monocyte_3_FCGR3A+
                                                     merged_df$rctd_Myeloid_c8_Monocyte_2_S100A9+
                                                     merged_df$rctd_Myeloid_c9_Macrophage_2_CXCL10+
                                                     merged_df$rctd_Cycling_Myeloid))
#Median myeloid calculation. 
Myeloid_stability<-median(merged_df$Myeloid_res)

##Normal epithelial resolution calculation
merged_df$Normal_epithelial_res<-abs(merged_df$rctd_Normal.Epithelial-(merged_df$rctd_Luminal.Progenitors+
                                                                       merged_df$rctd_Mature.Luminal+
                                                                       merged_df$rctd_Myoepithelial))
#Median epithelial calculation
Epithelial_stability<-median(merged_df$Normal_epithelial_res)

##Plasmablast resolution calculation
merged_df$Plasmablast_res<-abs(merged_df$rctd_Plasmablasts_major-merged_df$rctd_Plasmablasts_subset)
#Median plasmablast calculation
Plasmablast_stability<-median(merged_df$Plasmablast_res)

##PVL resolution calculation
merged_df$PVL_res<-abs(merged_df$rctd_PVL-(merged_df$rctd_PVL.Differentiated.s3+merged_df$rctd_PVL.Immature.s1+
                                             merged_df$rctd_PVL.Immature.s1+merged_df$rctd_PVL_Immature.s2+
                                             merged_df$rctd_Cycling.PVL))
#Median PVL calculaiton 
PVL_stability<-median(merged_df$PVL_res)

## T-cell resolution calculation
merged_df$T_cells_res<-abs(merged_df$rctd_T.cells-(merged_df$rctd_T_cells_c0_CD4._CCR7+
                                                   merged_df$rctd_T_cells_c10_NKT_cells_FCGR3A+
                                                   merged_df$rctd_T_cells_c11_MKI67+
                                                   merged_df$rctd_T_cells_c1_CD4._IL7R+
                                                   merged_df$rctd_T_cells_c2_CD4._T.regs_FOXP3+
                                                   merged_df$rctd_T_cells_c3_CD4._Tfh_CXCL13+
                                                   merged_df$rctd_T_cells_c4_CD8._ZFP36+
                                                   merged_df$rctd_T_cells_c5_CD8._GZMK+
                                                   merged_df$rctd_T_cells_c6_IFIT1+
                                                   merged_df$rctd_T_cells_c7_CD8._IFNG+
                                                   merged_df$rctd_T_cells_c8_CD8._LAG3+
                                                   merged_df$rctd_T_cells_c9_NK_cells_AREG))
#Median T-cell calculation
Tcell_stability<-median(merged_df$T_cells_res)
####PRODUCING FIGURES FOR INDIVIDUAL SPOTS and median differences for the given cell types
#Loading necessary libraries 
library(ggplot2)
library(tidyverse)
###MEDIAN DIFFERENCES FIGURE
##Making dataframe for the median differences for the given cell types (stability values)
Cell_types<-c("B-cells","CAFs","Cancer Epithelial","Endothelial","Epithelial",
              "Myeloid","Plasmablast","PVL","T-cells")
Stability_values<-c(Bcell_stability,CAF_stability,Cancer_Stability,
                    Endothelial_stability,Epithelial_stability,Myeloid_stability,
                    Plasmablast_stability,PVL_stability,Tcell_stability)
Stability_df<-data.frame(Cell_type=Cell_types,Stability=Stability_values)
#Plotting median differences for each cell type (stability values)
Stability<-Stability_df %>% 
  ggplot(aes(x=Cell_type,y=Stability)) +
  geom_col() +
  theme_bw() + 
  labs(title = "Stability in maintaining recall accuracy of RCTD at the subset phenotypic resolution for 
       TNBC",x="Major cell types",y="Stability metric") +
      theme(axis.text=element_text(size=12.5),
      axis.text.x = element_text(angle=90),
      axis.title.x = element_text(size = 15,face = "bold"),
      axis.title.y = element_text(size=15,face = "bold"),
      axis.text.x.bottom =element_text(vjust = 0.5),
      plot.title = element_text(size = 15,face = "bold"))
Stability
####FIGURE FOR DIFFERENCES IN EACH SPOT FOR MAJOR CELL TYPES
# Extract the last 9 columns
Spots<- merged_df[,64:72]
Spots_longer<-gather(Spots,key = "Cell type",value = "Difference")
Spot_plot<-Spots_longer %>% 
  ggplot(aes(x=`Cell type`,y=Difference)) +
  geom_boxplot() +
  theme_bw() + 
  labs(title = "Absolute difference between estimations by RCTD of major cell type and sum of its phenotypic subsets for each spot
       in context of TNBC",x="Major cell type",y="Absolute difference") +
  geom_jitter(color="black",size=0.4,alpha=0.9)
Spot_plot

Spot_plot + scale_x_discrete(labels = c("B-cells","CAFs", "Cancer Epithelial", "Endothelial", "Myeloid", "Normal Epithelial",
                                        "Plasmablast", "PVL", "T-cells"))


###WORKFLOW FOR OBTAINING RECALL CELL ACCURACY
##Loading in synthetic TNBC_major tissue pattern and obtaining relative composition
TNBC_major<-synthetic_visium_Wu_TNBC_major$relative_spot_composition
colnames(TNBC_major)[colnames(TNBC_major) == "name"] <-"spotID"
##Join TNBC_major and RCTD_major datasets by 'spotID'
Combined_df<-inner_join(TNBC_major,RCTD_major,by="spotID")
#Computing absolute difference of 'actual vs estimated'
#Cancer epithelial
Combined_df$Cancer_acc<-abs(Combined_df$Cancer.Epithelial-Combined_df$rctd_Cancer.Epithelial)
#Normal epithelial
Combined_df$Epithelial_acc<-abs(Combined_df$Normal.Epithelial-Combined_df$rctd_Normal.Epithelial)
#Endothelial
Combined_df$Endothelial_acc<-abs(Combined_df$Endothelial-Combined_df$rctd_Endothelial)
#PVL
Combined_df$PVL_acc<-abs(Combined_df$PVL-Combined_df$rctd_PVL)
#B-cells
Combined_df$Bcell_acc<-abs(Combined_df$B.cells-Combined_df$rctd_B.cells)
#CAFs
Combined_df$CAFs_acc<-abs(Combined_df$CAFs-Combined_df$rctd_CAFs)
#Myeloid
Combined_df$Myeloid_acc<-abs(Combined_df$Myeloid-Combined_df$rctd_Myeloid)
#T-cells
Combined_df$Tcell_acc<-abs(Combined_df$T.cells-Combined_df$rctd_T.cells)
#Plasmablasts
Combined_df$Plasmablast_acc<-abs(Combined_df$Plasmablasts-Combined_df$rctd_Plasmablasts)

#Boxplot showing accuracy for every spot
Accuracy<- Combined_df[,23:31]
Accuracy_longer<-gather(Accuracy,key="Cell type",value="Absolute Difference")
Accuracy_plot<-Accuracy_longer %>% 
  ggplot(aes(x=`Cell type`,y=`Absolute Difference`)) +
  geom_boxplot() +
  theme_bw() +
  labs(title="Absolute difference between actual major cell composition and predicted by RCTD for each spot within TNBC") +
  geom_jitter(color="black",size=0.4,alpha=0.9)

  Accuracy_plot + scale_x_discrete(labels=c("B-cells","CAFs","Cancer epithelial","Endothelial",
                                            "Myeloid","Normal Epithelial","Plasmablasts",
                                            "PVL","T-cells") )