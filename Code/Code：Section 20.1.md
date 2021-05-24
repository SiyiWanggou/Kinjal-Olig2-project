## **Code：Section 20.1**

### Section 20.1 Code branch for Kinjal's Olig2 project

```R
#------Senction 20------Olig2 Direction for Kinjal's Project
#Load packages
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(tidyverse)
library(dyno)
library(monocle)
library(ggpubr)
library(ggsignif)

setwd("/mnt/hgfs/share/Jerry_Project_Renewal/Ongoing/section20/")

#Load data
CGNP_P7_Math1_Cre_SmoM2_Sox2_reg <-
  readRDS(file = "/mnt/hgfs/share/Jerry_Project_Renewal/Saved_objects/section11_output/CC_pseudotime_Velocity_Projected_CGNP_P7_Math1_Cre_SmoM2_kcnb2_wt_Sox2_reg.Rds")

#Define cell cycle phase
CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$Phase_consistency <- 
  ifelse(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$S.Score_consistency1 < 0 & CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$G2M.Score_consistency1 < 0, "G1",
         ifelse(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$S.Score_consistency1 > CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$G2M.Score_consistency1, "S", "G2M"))

#------Neuronal differentiation Trajectory reconstruction P7-------
#loading Monocle object from Seruat object.
#Generate Sox2-MB-expression Matrix.
Mat_CC_Pseudotime_Sox2 <- CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@assays$RNA@counts
#Generate barcode file.
Barcode_CC_Pseudotime_Sox2 <- CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data
#Generate feature file.
Feature_CC_Pseudotime_Sox2 <- rownames(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg)
Feature_CC_Pseudotime_Sox2 <- as.data.frame(Feature_CC_Pseudotime_Sox2)
colnames(Feature_CC_Pseudotime_Sox2) <- c("gene_short_name")
rownames(Feature_CC_Pseudotime_Sox2) <- Feature_CC_Pseudotime_Sox2$gene_short_name
#Establish the monocle2 object
pd<-new("AnnotatedDataFrame",data = Barcode_CC_Pseudotime_Sox2)
fd<-new("AnnotatedDataFrame",data = Feature_CC_Pseudotime_Sox2)
Monocle_CC_Pseudotime_Sox2 <- newCellDataSet(Mat_CC_Pseudotime_Sox2,
                                             phenoData = pd,
                                             featureData = fd,
                                             lowerDetectionLimit = 0.1,
                                             expressionFamily = negbinomial.size()) 
#Estimate Size Factors and Dispersions.
Monocle_CC_Pseudotime_Sox2 <- estimateSizeFactors(Monocle_CC_Pseudotime_Sox2)
Monocle_CC_Pseudotime_Sox2 <- estimateDispersions(Monocle_CC_Pseudotime_Sox2)
Monocle_CC_Pseudotime_Sox2 <- detectGenes(Monocle_CC_Pseudotime_Sox2,min_expr = 0.1) 
expressed_genes <- row.names(subset(fData(Monocle_CC_Pseudotime_Sox2),num_cells_expressed >=10)) 
Monocle_CC_Pseudotime_Sox2 <- setOrderingFilter(Monocle_CC_Pseudotime_Sox2,expressed_genes)

#Generate Neurogeneisis associated genes for trajectory reconstruction.
GO_neurogeneisis <- read.delim("/mnt/hgfs/share/Jerry_Project_Renewal/Saved_objects/prepared_objects/GO_term_20200211_161335.txt",header = TRUE)
Gene_list <- unique(GO_neurogeneisis$Symbol)

#Pre-identify DE genes in Sox2-MB-cells according to seurat_clusters.
clustering_DEG_genes <- differentialGeneTest(Monocle_CC_Pseudotime_Sox2[expressed_genes,],
                                             fullModelFormulaStr = '~ projected_sample_age + Phase_consistency',
                                             reducedModelFormulaStr = '~ Phase_consistency',
                                             cores = 12, verbose = TRUE)

#Function for model selection
Pseudotime_Trajectory_Reconstruction <- function(x,y,g1,g2){
  HSMM_data <- x
  Sig_DE_genes <- subset(clustering_DEG_genes,qval < y)
  Gene_index <- row.names(Sig_DE_genes[which(rownames(Sig_DE_genes) %in% Gene_list),])
  HSMM_data <- setOrderingFilter(HSMM_data,Gene_index)
  HSMM_data <- reduceDimension(HSMM_data, max_components = 2, num_dim = 6, 
                               norm_method = "log", pseudo_expr = 1,
                               reduction_method = 'DDRTree', verbose = T, 
                               residualModelFormulaStr = "~ Phase_consistency + num_genes_expressed",cores = detectCores() - 2)
  HSMM_data <- orderCells(HSMM_data)
  #Plot Trajectory
  plot_a <- plot_cell_trajectory(HSMM_data, max_components = 2, color_by = "projected_sample_age", show_tree = FALSE,
                                 cell_size = 1,cell_link_size = 1, state_number_size = 0.1,show_branch_points = TRUE) + 
    scale_colour_manual(values = c("#FF0000","#FF8000","#FFFF00","#80FF00","#00FF00","00FF80","#00FFFF","#0080FF","#0000FF","#7F00FF","#FF00FF","#FF007F","#808080")) +
    facet_wrap(~sample) +
    theme(legend.text=element_text(size=3)) + 
    theme(legend.position="top",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  #plot branch1 markers
  target<-rownames(fData(HSMM_data)[which(fData(HSMM_data)$gene_short_name==g1),])
  pData(HSMM_data)$g1 <- exprs(HSMM_data)[target,]
  plot_b <- plot_cell_trajectory(HSMM_data, color_by = "g1", 
                                 cell_size = 1,cell_link_size = 1, state_number_size = 0.1,show_branch_points = FALSE) + 
    scale_colour_gradient2(low="white",mid="Magenta2",high="Red",midpoint = 10,limits = c(0,20)) +
    facet_wrap(~sample) +
    theme(legend.position="top",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  #plot branch2 markers  
  target<-rownames(fData(HSMM_data)[which(fData(HSMM_data)$gene_short_name==g2),])
  pData(HSMM_data)$g2 <-exprs(HSMM_data)[target,]
  plot_c <- plot_cell_trajectory(HSMM_data, color_by = "g2", 
                                 cell_size = 1,cell_link_size = 1, state_number_size = 0.1,show_branch_points = FALSE) + 
    scale_colour_gradient2(low="white",mid="Magenta2",high="Red",midpoint = 10,limits = c(0,20)) +
    facet_wrap(~sample) +
    theme(legend.position="top",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  #Merge the plots
  plot_merge <- ggarrange(plot_a,plot_b,plot_c,nrow = 3)
  results <- list(HSMM_data,plot_merge) #Monocle object stored in [[1]], plots stored in [[2]].
  return(results)
}

#Generate the multiple test results from qval < 0.1, 0.05, 0.01, 0.005, 0.001.
Monocle_CC_Pseudotime_Sox2_01 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.1,"Stmn2","Olig2")
Monocle_CC_Pseudotime_Sox2_0075 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.075,"Stmn2","Olig2")
Monocle_CC_Pseudotime_Sox2_005 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.05,"Stmn2","Olig2")
Monocle_CC_Pseudotime_Sox2_0025 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.025,"Stmn2","Olig2")
Monocle_CC_Pseudotime_Sox2_001 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.01,"Stmn2","Olig2")
Monocle_CC_Pseudotime_Sox2_0005 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.005,"Stmn2","Olig2")
Monocle_CC_Pseudotime_Sox2_0001 <- Pseudotime_Trajectory_Reconstruction(Monocle_CC_Pseudotime_Sox2,0.001,"Stmn2","Olig2")

#Model selection
Cairo(file="Model_Seletion_CC_Pseudotime_Sox2.png",type="png",units="in",bg="white",width=24,height=12,pointsize=5,dpi=300)
ggarrange(Monocle_CC_Pseudotime_Sox2_01[[2]],Monocle_CC_Pseudotime_Sox2_0075[[2]],Monocle_CC_Pseudotime_Sox2_005[[2]],
          Monocle_CC_Pseudotime_Sox2_0025[[2]],Monocle_CC_Pseudotime_Sox2_001[[2]],Monocle_CC_Pseudotime_Sox2_0005[[2]],
          Monocle_CC_Pseudotime_Sox2_0001[[2]],
          labels = c("q < 0.1","q < 0.075","q < 0.05","q < 0.025","q < 0.01","q < 0.005","q < 0.001"),
          ncol = 7)
dev.off()
```

![Model_Seletion_CC_Pseudotime_Sox2.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Model_Seletion_CC_Pseudotime_Sox2.png?raw=true)

```R
#Accept "q < 0.01" for model reconstruction
Monocle_CC_Pseudotime_Sox2 <- Monocle_CC_Pseudotime_Sox2_001[[1]]

#Identify root branch
Monocle_CC_Pseudotime_Sox2 <- orderCells(Monocle_CC_Pseudotime_Sox2, root_state = 4)

#plot the pseudotime
Average <- (max(pData(Monocle_CC_Pseudotime_Sox2)$Pseudotime)+min(pData(Monocle_CC_Pseudotime_Sox2)$Pseudotime))/2
Cairo(file="Trajectory_Monocle_CC_Pseudotime_Sox2_pseudotime.png",type="png",units="in",bg="white",width=4,height=4.5,pointsize=10,dpi=300)
plot_cell_trajectory(Monocle_CC_Pseudotime_Sox2, 1, 2, color_by = "Pseudotime", show_tree = TRUE,
                     cell_size = 1.5,cell_link_size = 0.5, state_number_size = 0.1,show_branch_points = FALSE) + 
  scale_colour_gradient2(low="#00FF00",mid="#FF8000",high="#FF0000",midpoint = Average) +
  facet_wrap(~sample) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme(legend.position="top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

#Identify cell fate directions
pData(Monocle_CC_Pseudotime_Sox2)$cellfate <- dplyr::recode(pData(Monocle_CC_Pseudotime_Sox2)$State,
                                                            "4" = "Pre_branch",
                                                            "1" = "Olig2_direction", "2" = "Olig2_direction", "3" = "Olig2_direction",
                                                            "5" = "Neurod1_direction", "6" = "Neurod1_direction", "7" = "Neurod1_direction")


#Plot RNA velocity onto Monocle2 DDRTree trajectory
embeddings <- t(Monocle_CC_Pseudotime_Sox2@reducedDimS)
colnames(embeddings) <- c("Dim1","Dim2")
pData(Monocle_CC_Pseudotime_Sox2)
CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data <- cbind(pData(Monocle_CC_Pseudotime_Sox2),embeddings)

#Prepare embedding neurodifferentiation
Emb <- as.matrix(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data[,c("Dim1","Dim2")])

#Assign the cell colors
cell_colors <- cut(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$Pseudotime,
                   quantile(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$Pseudotime,probs = seq(0,1,0.1)),
                   labels = viridis::plasma(10))
names(cell_colors) <- rownames(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data)

#Total Arrow of RNA velocity
Cairo(file="Velocity_arrow_on_cells.png",type="png",units="in",bg="white",width=5,height=5,pointsize=5,dpi=300)
show.velocity.on.embedding.cor(emb =Emb,vel = Tool(object = CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),n = 200,
                               scale = "sqrt",cell.colors = ac(x = cell_colors, alpha = 1.0),
                               cex = 4, arrow.scale = 1.5, show.grid.flow = FALSE, min.grid.cell.mass = 1, grid.n = 50,
                               arrow.lwd = 1.25,do.par = FALSE, cell.border.alpha = 0,n.cores = 14)
dev.off()
Cairo(file="Velocity_field_on_cells.png",type="png",units="in",bg="white",width=5,height=5,pointsize=5,dpi=300)
show.velocity.on.embedding.cor(emb =Emb,vel = Tool(object = CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),n = 200,
                               scale = "sqrt",cell.colors = ac(x = cell_colors, alpha = 1.0),
                               cex = 4, arrow.scale = 1.5, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 50,
                               arrow.lwd = 1.25,do.par = FALSE, cell.border.alpha = 0,n.cores = 14)
dev.off()
```

![Velocity_arrow_on_cells.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_arrow_on_cells.png?raw=true)

![Velocity_field_on_cells.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_field_on_cells.png?raw=true)



```R
#Plot maker gene velocity
#Neurod1
Cairo(file="Velocity_Markers_Neurod1.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Neurod1", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Olig2
Cairo(file="Velocity_Markers_Olig2.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Olig2", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Pdgfra
Cairo(file="Velocity_Markers_Pdgfra.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Pdgfra", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Sox10
Cairo(file="Velocity_Markers_Sox10.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Sox10", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Mki67
Cairo(file="Velocity_Markers_Mki67.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Mki67", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Top2a
Cairo(file="Velocity_Markers_Top2a.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Top2a", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Cdk1
Cairo(file="Velocity_Markers_Cdk1.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = Emb, #cell.colors = cell_colors_GCP,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Cdk1", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
```

![Velocity_Markers_Olig2.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Olig2.png?raw=true)

![Velocity_Markers_Neurod1.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Neurod1.png?raw=true)

![Velocity_Markers_Mki67.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Mki67.png?raw=true)

![Velocity_Markers_Top2a.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Top2a.png?raw=true)

![Velocity_Markers_Cdk1.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Cdk1.png?raw=true)

![Velocity_Markers_Pdgfra.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Pdgfra.png?raw=true)

![Velocity_Markers_Sox10.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Sox10.png?raw=true)



```R
#------Cell cycle Trajectory reconstruction P7-------
#Reconstruction of cell cycle trajectory at P7
DefaultAssay(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg) <- "RNA"
Dyno_P7 <- wrap_expression(counts = t(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@assays$RNA@counts),
                           expression = t(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@assays$RNA@data))
guidelines <- guidelines_shiny(Dyno_P7)
addTaskCallback(function(...) {set.seed(123);TRUE})
P7_model_slice <- infer_trajectories(Dyno_P7,c("slingshot"),seed = 123,verbose = TRUE)
P7_model_slice_results <- P7_model_slice[[6]][[1]]
dimred <- CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data[,c("S.Score_consistency1","G2M.Score_consistency1")]
#Results for slice model
#Plot Trajectory
P7_model_slice_results <- P7_model_slice_results %>% add_root(root_milestone_id = "4") #Select the root for pseudotime
simplified <- simplify_trajectory(P7_model_slice_results)
overall_feature_importances <- dynfeature::calculate_overall_feature_importance(P7_model_slice_results,expression_source = Dyno_P7$expression) #Generate feature genes along trajectory overall
features <- overall_feature_importances %>% top_n(40, importance) %>% pull(feature_id) #Generate feature genes along trajectory overall
Cairo(file="Dyno_P7_slingshot_Trajectory.png",type="png",units="in",bg="white",width=16,height=6,pointsize=10,dpi=300)
patchwork::wrap_plots(
  plot_dimred(P7_model_slice_results,expression_source = Dyno_P7$expression,grouping = group_onto_nearest_milestones(P7_model_slice_results),dimred = dimred) + ggtitle("Milestones") +
    coord_fixed(ratio = 0.7), #Trajectory marked with milestones.
  plot_dimred(P7_model_slice_results,'pseudotime',pseudotime = calculate_pseudotime(P7_model_slice_results),dimred = dimred) + ggtitle("Pseudotime") +
    scale_colour_gradient2(low="purple",mid="#FF007F",high="Yellow",midpoint = 10) + ggtitle("Pseudotime") +
    coord_fixed(ratio = 0.7), #Trajectory marked with pseudotime
  plot_dimred(P7_model_slice_results,expression_source = Dyno_P7$counts, dimred = dimred, feature_oi = "Ung") + coord_fixed(ratio = 0.7) +
    scale_colour_gradient2(low="white",mid="Magenta2",high="red",midpoint = 5) + ggtitle("Ung"), #Trajectory marked with Tubb3 expression.
  plot_dimred(P7_model_slice_results,expression_source = Dyno_P7$counts, dimred = dimred, feature_oi = "Pcna") + coord_fixed(ratio = 0.7) +
    scale_colour_gradient2(low="white",mid="Magenta2",high="red",midpoint = 5) + ggtitle("Pcna"), #Trajectory marked with Tubb3 expression.
  plot_dimred(P7_model_slice_results,expression_source = Dyno_P7$counts, dimred = dimred, feature_oi = "Mki67") + coord_fixed(ratio = 0.7) +
    scale_colour_gradient2(low="white",mid="Magenta2",high="red",midpoint = 50) + ggtitle("Mki67"), 
  plot_dimred(P7_model_slice_results,expression_source = Dyno_P7$counts, dimred = dimred, feature_oi = "Top2a") + coord_fixed(ratio = 0.7) +
    scale_colour_gradient2(low="white",mid="Magenta2",high="red",midpoint = 41) + ggtitle("Top2a"), # Trajectory marked with Sox2 expression.
  #Trajectory marked with Mki67 expression.
  ncol = 6
)
dev.off()
```

![Dyno_P7_slingshot_Trajectory.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Dyno_P7_slingshot_Trajectory.png?raw=true)

```R
#Generate the pseudotime along cell cycle at P7
P7_pseudotime_cell_cycle_slice <- calculate_pseudotime(P7_model_slice_results)
CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$Pseudotime_cell_cycle_regain <- P7_pseudotime_cell_cycle_slice

#Plot RNA velocity of P7 on cell cycling score embedding
#Assign the cell colors
cell_colors_cellcycle <- cut(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$Pseudotime_cell_cycle,
                             quantile(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data$Pseudotime_cell_cycle,probs = seq(0,1,0.12)),
                             labels = viridis::plasma(8))
names(cell_colors_cellcycle) <- rownames(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data)

#Total Arrow of RNA velocity in cell cycling based on slingshot pseudotime
dimred <- as.matrix(dimred)
Cairo(file="P7_Velocity_Cellcycle_on_Cell_Cycle_pseudotime_arrow_on_cells.png",type="png",units="in",bg="white",width=5,height=5,pointsize=5,dpi=300)
show.velocity.on.embedding.cor(emb = dimred,vel = Tool(object = CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, 
                                                       slot = "RunVelocity"),n = 200,
                               scale = "sqrt",cell.colors = ac(x = cell_colors_cellcycle, alpha = 1.0),
                               cex = 3.0, arrow.scale = 0.5, show.grid.flow = FALSE, min.grid.cell.mass = 1, grid.n = 50,
                               arrow.lwd = 1.0,do.par = FALSE, cell.border.alpha = 0,n.cores = 14)
dev.off()
Cairo(file="P7_Velocity_Cellcycle_on_Cell_Cycle_pseudotime_arrow_on_fields.png",type="png",units="in",bg="white",width=5,height=5,pointsize=5,dpi=300)
show.velocity.on.embedding.cor(emb = dimred,vel = Tool(object = CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, 
                                                       slot = "RunVelocity"),n = 200,
                               scale = "sqrt",cell.colors = ac(x = cell_colors_cellcycle, alpha = 1.0),
                               cex = 3.0, arrow.scale = 0.5, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 50,
                               arrow.lwd = 1.0,do.par = FALSE, cell.border.alpha = 0,n.cores = 14)
dev.off()
```

![P7_Velocity_Cellcycle_on_Cell_Cycle_pseudotime_arrow_on_cells.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/P7_Velocity_Cellcycle_on_Cell_Cycle_pseudotime_arrow_on_cells.png?raw=true)

![P7_Velocity_Cellcycle_on_Cell_Cycle_pseudotime_arrow_on_fields.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/P7_Velocity_Cellcycle_on_Cell_Cycle_pseudotime_arrow_on_fields.png?raw=true)

```R
#Marker gene RNA velocity on cell cycle embedding.
#Cdk1
Cairo(file="Velocity_Markers_Cdk1_CC.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = dimred,
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Cdk1", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Pcna
Cairo(file="Velocity_Markers_Pcna_CC.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = dimred, 
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Pcna", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Mki67
Cairo(file="Velocity_Markers_Mki67_CC.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = dimred, 
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Mki67", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Neurod1
Cairo(file="Velocity_Markers_Neurod1_CC.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = dimred, 
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Neurod1", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
#Olig2
Cairo(file="Velocity_Markers_Olig2_CC.png",type="png",units="in",bg="white",width=22,height=6,pointsize=25,dpi=300)
gene.relative.velocity.estimates(GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "spliced"),
                                 GetAssayData(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg,slot = "data", assay = "unspliced"),
                                 deltaT = 1, kCells = 7, fit.quantile = 0.02, n.cores = 14,
                                 cell.emb = dimred, 
                                 old.fit = Tool(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg, slot = "RunVelocity"),diagonal.quantiles = FALSE,
                                 show.gene = "Olig2", expression.gradient = RColorBrewer::brewer.pal(9,"Greens"),residual.gradient = NULL)
dev.off()
```

![Velocity_Markers_Mki67_CC.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Mki67_CC.png?raw=true)

![Velocity_Markers_Cdk1_CC.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Cdk1_CC.png?raw=true)

![Velocity_Markers_Pcna_CC.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Pcna_CC.png?raw=true)

![Velocity_Markers_Olig2_CC.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Olig2_CC.png?raw=true)

![Velocity_Markers_Neurod1_CC.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Velocity_Markers_Neurod1_CC.png?raw=true)

```R
#Plot cell cycle phases against the two branches.
#loading marker genes
P7_metadata <- cbind(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data,
                     t(as.matrix(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@assays$RNA@data[c("Olig2","Mki67","Top2a","Cdk1","Neurod1"),])))
Olig2_direction_metadata <- subset(P7_metadata,cellfate == "Olig2_direction" | cellfate == "Pre_branch")
Neurod1_direction_metadata <- subset(P7_metadata,cellfate == "Neurod1_direction" | cellfate == "Pre_branch")
#Recombinding the dataset
Olig2_direction_metadata$cellfate_1 <- c(rep("Olig2",nrow(Olig2_direction_metadata)))
Neurod1_direction_metadata$cellfate_1 <- c(rep("Neurod1",nrow(Neurod1_direction_metadata)))
merge_metadata <- rbind(Olig2_direction_metadata,Neurod1_direction_metadata)
#Loess fitting.
q_Combined <- ggplot(data=merge_metadata,aes(x=Pseudotime,y=Pseudotime_cell_cycle_regain,color = cellfate_1)) +
  stat_smooth(method = 'loess',level = 0.75, span = 0.2,size = 1.5) + 
  coord_cartesian(xlim=c(0,7)) +
  scale_color_manual(values = c("Olig2" = "Red", "Neurod1" = "black")) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Neuronal differentiation pseudotime",y="Cell cycle pseudotime")
#Olig2 expression along the two trajectories
q_Olig2 <- ggplot(data=merge_metadata,aes(x=Pseudotime,y=Olig2,color = cellfate_1)) +
  stat_smooth(method = 'loess',level = 0.3, span = 0.2,size = 1.5) + 
  scale_color_manual(values = c("Olig2" = "Red", "Neurod1" = "black")) +
  coord_cartesian(xlim=c(0,7)) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Neuronal differentiation pseudotime",y="Olig2 expression")
#Neurod1 expression along the two trajectories
q_Neurod1 <- ggplot(data=merge_metadata,aes(x=Pseudotime,y=Neurod1,color = cellfate_1)) +
  stat_smooth(method = 'loess',level = 0.3, span = 0.2,size = 1.5) + 
  scale_color_manual(values = c("Olig2" = "Red", "Neurod1" = "black")) +
  coord_cartesian(xlim=c(0,7)) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Neuronal differentiation pseudotime",y="Neurod1 expression")

#Cell cycling gene expression along the two trajectories
q_Mki67 <- ggplot(data=merge_metadata,aes(x=Pseudotime,y=Mki67,color = cellfate_1)) +
  stat_smooth(method = 'loess',level = 0.3, span = 0.2,size = 1.5) + 
  scale_color_manual(values = c("Olig2" = "Red", "Neurod1" = "black")) +
  coord_cartesian(xlim=c(0,7)) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Neuronal differentiation pseudotime",y="Mki67 pseudotime")
q_Top2a <- ggplot(data=merge_metadata,aes(x=Pseudotime,y=Top2a,color = cellfate_1)) +
  stat_smooth(method = 'loess',level = 0.3, span = 0.2,size = 1.5) + 
  scale_color_manual(values = c("Olig2" = "Red", "Neurod1" = "black")) +
  coord_cartesian(xlim=c(0,7)) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Neuronal differentiation pseudotime",y="Top2a expression")
q_Cdk1 <- ggplot(data=merge_metadata,aes(x=Pseudotime,y=Cdk1,color = cellfate_1)) +
  stat_smooth(method = 'loess',level = 0.3, span = 0.2,size = 1.5) + 
  scale_color_manual(values = c("Olig2" = "Red", "Neurod1" = "black")) +
  theme_bw() +
  theme(legend.position="top",text=element_text(size=25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = "white")) +
  labs(x="Neuronal differentiation pseudotime",y="Cdk1 Expression")

Cairo(file="Olig2_Neurod1_Expression.png",type="png",units="in",bg="white",width=8,height=8,pointsize=5,dpi=300)
patchwork::wrap_plots(q_Olig2, q_Neurod1,ncol = 1)
dev.off()

Cairo(file="Pseudotime_combination_Mki67_Top2a_Cdk1_Expression.png",type="png",units="in",bg="white",width=8,height=16,pointsize=5,dpi=300)
patchwork::wrap_plots(q_Combined, q_Mki67, q_Top2a, q_Cdk1, ncol = 1)
dev.off()
```

![Olig2_Neurod1_Expression.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Olig2_Neurod1_Expression.png?raw=true)

![Pseudotime_combination_Mki67_Top2a_Cdk1_Expression.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Pseudotime_combination_Mki67_Top2a_Cdk1_Expression.png?raw=true)

```R
#Identify Olig2 confirmed and Neurod1 confirmed Sox2 cells
#Comparison between cell fate confirmed Sox2 cells
P7_metadata <- cbind(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@meta.data,
                     t(as.matrix(CGNP_P7_Math1_Cre_SmoM2_Sox2_reg@assays$RNA@data[c("Olig2","Mki67","Top2a","Cdk1","Neurod1"),])))
P7_metadata$Cell_fate_confirmation <- ifelse(P7_metadata$Pseudotime >= 5.2 & P7_metadata$cellfate == "Olig2_direction", 
                                             "Cellfate_Confirmed_Olig2",ifelse(P7_metadata$Pseudotime >= 4.9 & P7_metadata$cellfate == "Neurod1_direction", 
                                                                               "Cellfate_Confirmed_Neurod1","Others"))
#Mki67_comparison
Cairo(file="Mki67_comparison.png",type="png",units="in",bg="white",width=4,height=2,pointsize=10,dpi=300)
ggplot(subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"),
       aes(x=Cell_fate_confirmation,y=Mki67, fill = Cell_fate_confirmation)) +
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Mki67") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#Top2a_comparison
Cairo(file="Top2a_comparison.png",type="png",units="in",bg="white",width=4,height=2,pointsize=10,dpi=300)
ggplot(subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"),
       aes(x=Cell_fate_confirmation,y=Top2a, fill = Cell_fate_confirmation)) + 
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Top2a") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#Cdk1_comparison
Cairo(file="Cdk1_comparison.png",type="png",units="in",bg="white",width=4,height=2,pointsize=10,dpi=300)
ggplot(subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"),
       aes(x=Cell_fate_confirmation,y=Cdk1, fill = Cell_fate_confirmation)) +
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Cdk1") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#Statstics for marker gene comparsion
cells_id <- rownames(subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"))
marker_genes <- c("Mki67","Top2a","Cdk1")
Monocle_Cellfate_Confirmed_Sox2 <- Monocle_CC_Pseudotime_Sox2[marker_genes,cells_id]
pData(Monocle_Cellfate_Confirmed_Sox2)$cell_fate_confirmation <- subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1")$Cell_fate_confirmation
clustering_DEG_genes_cell_fate_confirmed <- differentialGeneTest(Monocle_Cellfate_Confirmed_Sox2,
                                             fullModelFormulaStr = '~ cell_fate_confirmation',
                                             cores = 12, verbose = TRUE)
clustering_DEG_genes_cell_fate_confirmed
#Cell cycle pseudotime comparison
Cairo(file="Pseudotime_cell_cycle_regain_comparison.png",type="png",units="in",bg="white",width=4,height=2,pointsize=10,dpi=300)
ggplot(subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"),
       aes(x=Cell_fate_confirmation,y=Pseudotime_cell_cycle_regain, fill = Cell_fate_confirmation)) + 
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Pseudotime_cell_cycle_regain") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
a <- aov(Pseudotime_cell_cycle_regain~Cell_fate_confirmation,subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"))
b <- wilcox.test(Pseudotime_cell_cycle_regain~Cell_fate_confirmation,subset(P7_metadata,Cell_fate_confirmation == "Cellfate_Confirmed_Olig2" | Cell_fate_confirmation == "Cellfate_Confirmed_Neurod1"))
b
summary(a)
```

![Pseudotime_cell_cycle_regain_comparison.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Pseudotime_cell_cycle_regain_comparison.png?raw=true)

![Mki67_comparison.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Mki67_comparison.png?raw=true)

![Top2a_comparison.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Top2a_comparison.png?raw=true)

![Cdk1_comparison.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Cdk1_comparison.png?raw=true)

```R
#Comparison between Olig2 cell fate confirmed Sox2 cells and other Sox2 cells
P7_metadata$Cell_fate_confirmation <- ifelse(P7_metadata$Pseudotime >= 5.2 & P7_metadata$cellfate == "Olig2_direction", 
                                             "Cellfate_Confirmed_Olig2","Others")
#Mki67_comparison
Cairo(file="Mki67_comparison_with_others.png",type="png",units="in",bg="white",width=4,height=2.5,pointsize=10,dpi=300)
ggplot(P7_metadata,
       aes(x=Cell_fate_confirmation,y=Mki67, fill = Cell_fate_confirmation)) + 
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Mki67") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#Top2a_comparison
Cairo(file="Top2a_comparison_with_others.png",type="png",units="in",bg="white",width=4,height=2.5,pointsize=10,dpi=300)
ggplot(P7_metadata,
       aes(x=Cell_fate_confirmation,y=Top2a, fill = Cell_fate_confirmation)) + 
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Top2a") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#Cdk1_comparison
Cairo(file="Cdk1_comparison_with_others.png",type="png",units="in",bg="white",width=4,height=2.5,pointsize=10,dpi=300)
ggplot(P7_metadata,
       aes(x=Cell_fate_confirmation,y=Cdk1, fill = Cell_fate_confirmation)) + 
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Cdk1") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
#Statstics for marker gene comparsion
marker_genes <- c("Mki67","Top2a","Cdk1")
Monocle_Cellfate_Confirmed_Sox2 <- Monocle_CC_Pseudotime_Sox2[marker_genes,]
pData(Monocle_Cellfate_Confirmed_Sox2)$cell_fate_confirmation <- P7_metadata$Cell_fate_confirmation
clustering_DEG_genes_cell_fate_confirmed <- differentialGeneTest(Monocle_Cellfate_Confirmed_Sox2,
                                                                 fullModelFormulaStr = '~ cell_fate_confirmation',
                                                                 cores = 12, verbose = TRUE)
clustering_DEG_genes_cell_fate_confirmed
#Cell cycle pseudotime comparison
Cairo(file="Pseudotime_cell_cycle_regain_comparison_with_Others.png",type="png",units="in",bg="white",width=4,height=2.5,pointsize=10,dpi=300)
ggplot(P7_metadata,
       aes(x=Cell_fate_confirmation,y=Pseudotime_cell_cycle_regain, fill = Cell_fate_confirmation)) + 
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width = 0.1,fill="white")+
  labs(title="Violin Plot",x="Genotype",y="Pseudotime_cell_cycle_regain") +
  scale_fill_brewer(palette = "RdBu") +
  theme_classic()
dev.off()
a <- aov(Pseudotime_cell_cycle_regain~Cell_fate_confirmation,P7_metadata)
b <- wilcox.test(Pseudotime_cell_cycle_regain~Cell_fate_confirmation,P7_metadata)
b
summary(a)
```

![Pseudotime_cell_cycle_regain_comparison_with_Others.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Pseudotime_cell_cycle_regain_comparison_with_Others.png?raw=true)

![Mki67_comparison_with_others.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Mki67_comparison_with_others.png?raw=true)

![Top2a_comparison_with_others.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Top2a_comparison_with_others.png?raw=true)

![Cdk1_comparison_with_others.png](https://github.com/SiyiWanggou/Kinjal-Olig2-project/blob/main/Results/Cdk1_comparison_with_others.png?raw=true)