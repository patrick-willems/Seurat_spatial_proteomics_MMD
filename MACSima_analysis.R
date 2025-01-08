## Script name: MACSima analysis
##
## Use the Seurat to analyse MACSima fluorescence microscopy data of a Moyamoya biopsy sample.
##
## Author: Dr. Patrick Willems
##
## Email: patwille.willems@UGent.be



#Load all libraries and set the working directory.
library(Seurat)
library(dplyr)
library(ggsci)
library(scCustomize)
library(ggplot2)
library(patchwork)
library(tidyr)
setwd("C:/Users/Installer/Documents/DATA/MACSima/ONLINE/")

# Processed MASCima fluorescence intensities were re-formatted in a 10X-like object with
# features = channels, barcodes = segmented cells and matrix = average fluorescence intensities
data <- Read10X("10X_input/",gene.column = 1)

# Load metadata from the spatial proteomics image analysis
metaData <- read.table('metadata.txt',sep='\t',row.names = 1,header=T) #Metadata

#Generate a Seurat object with metadata
seuratObj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 2) %>%
  AddMetaData(metaData)

# Remove cells presenting outlier intensity values for certain channels and cells
# presenting acquisition bleaching or residing in out-of-focus and tissue folding areas. 
filteredObj <- subset(x=seuratObj, subset = artifact_annotation == "Not_assigned" & is.na(seuratObj@meta.data$outlier))

# This goes from 3,322 total cells to 2,362 high-quality cells.
dim(seuratObj)
dim(filteredObj)

# Follow the CITE-seq normalisation procedure and standard Seurat processing
filteredObj <- NormalizeData(filteredObj,normalization.method='CLR',margin = 2) %>%
  FindVariableFeatures() %>%
  ScaleData()  %>%
  RunPCA(npcs=50)

# Determine number of components to be used.
ElbowPlot(filteredObj,ndims = 50)
DimHeatmap(filteredObj,dims = 1:12,balanced = TRUE,cells=500)
DimHeatmap(filteredObj,dims = 13:24,balanced = TRUE,cells=500)

# Set to 12 components.
filteredObj <- FindNeighbors(filteredObj,dims=1:12) %>%
  FindClusters(resolution=0.2,graph.name="RNA_snn") %>%
  RunUMAP(dims=1:12,label=T,min.dist = 0.1) 

# UMAP plot.
my_colors = c('#f0e443','#d56014','#069f73','#0173b4','#cc79a7')
DimPlot(filteredObj,label=T,label.size=5,pt.size=.7,reduction = "umap",cols = my_colors) + NoLegend() + NoAxes()

# Color by segmentation (immunecell, red blood cell or other)
DimPlot(filteredObj, label=F, pt.size=.5, group.by="original_segmentation_mask", reduction = "umap")

# Color by pre-defined region (media, intima, adventitia, ..)
DimPlot(filteredObj, label=F, pt.size=.5, group.by="region_annotation", reduction = "umap")

# Plot a feature.
palette <- c("#808080","cornsilk3", "#FFFF00","red3")
FeaturePlot_scCustom(filteredObj,na_color="grey",order=T,colors=palette,pt.size=0.8,features = "Desmin-C-REA1134",label=F,reduction = "umap")

# Retrieve most disriminative markers per cluster.
markers <- FindAllMarkers(filteredObj, test.use = "roc", only.pos = TRUE)
#write.table(markers,'allmarkers_roc.txt',sep='\t')

# Get the top 5 markers (columns) per cluster (rows)
top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(power), desc(abs(avg_log2FC))) %>%
  slice_head(n = 5) %>%
  arrange(cluster, desc(power))

# Make combined violin plots: 5 rows = clusters and 5 genes = columns
plots <- list()
for (i in 1:nrow(top_markers)) {
  gene <- top_markers$gene[i]
  cluster <- top_markers$cluster[i]
  adj_p_val <-formatC(top_markers$power[i], digits = 3)
  plot <- VlnPlot(filteredObj, features = gene, pt.size = 0.1) +
    labs(title = paste(gene, ", ",cluster,":",adj_p_val,sep='')) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank()   # Remove y-axis title
    )
  
  # Add the plot to the list
  plots[[paste0("Cluster_", cluster, "_Gene_", gene)]] <- plot
}
wrap_plots(plots, ncol = 5)

# Focus on segmented immunecells (N=223)
immuneObj <- subset(x = filteredObj, subset = original_segmentation_mask == "segmentation_mask_immune")
dim(immuneObj)

# Pie charts according region_annotation
region_colors = c('#069f73','#f0e443','#d56014','#0173b4','firebrick','grey')
immune_per_region <-  as.data.frame(table(immuneObj$region_annotation))
colnames(immune_per_region) <- c("Region", "Count")
immune_per_region$Region <- factor(immune_per_region$Region, levels = immune_per_region$Region[order(-immune_per_region$Count)])
ggplot(immune_per_region, aes(x = "", y = Count, fill = factor(Region))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = region_colors) +
  labs(title = "Gated Immune Cells", fill = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())


# Plot the markers per region
region_colors <- c('#f0e443','#069f73','#d56014')
names(region_colors) <- c('intima', 'tunica_adventitia', 'lumen')
immune_markers <- c("Syk-C-REA111","HLADR-C-REAL550","CD8a-C-REA1024","CD66b-C-REA306","CD57-C-REA769","CD45RB-C-REA11","CD45RA-C-REAL164","CD45-C-5B1","CD43-C-REA833","CD2-C-REA1130","CD162-C-REA319","CD15-C-VIMC6")

# Fetch data and put in a long format for ggplot2-based plotting
immune_data <- FetchData(immuneObj, vars = immune_markers)
immune_data$Region <- immuneObj$region_annotation
immune_data_long <- immune_data %>%
  pivot_longer(cols = all_of(immune_markers), names_to = "Marker", values_to = "Expression") %>%
  filter(Region != "Not_assigned")

for(marker in immune_markers){
  print(marker)
  # Filter the data for the current marker
  immune_data_marker <- immune_data_long %>%
    filter(Marker == marker) %>%
    filter(Region %in% c('tunica_adventitia', 'intima', 'lumen'))  # Filter for the regions of interest
  
  # Plot for the current marker across all regions
  ggplot(immune_data_marker, aes(x = Marker, y = Expression, fill = Region)) +
    geom_violin(aes(fill = Region), trim = FALSE, alpha = 1, scale = "width", position = position_dodge(width = 1)) +  # Violin plot with dodge
    geom_boxplot(aes(group = Region), width = 0.2, color = "black", fill = "grey", outlier.shape = NA, coef = 0, position = position_dodge(width = 1)) +  # Boxplot with dodge
    geom_jitter(aes(group = Region), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.1)) +  # Jitter with dodge width
    scale_fill_manual(values = region_colors) +  # Use the color palette for the marker
    theme_minimal() +
    labs(title = paste("Marker:", marker), x = "Immune Marker", y = "Expression") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.background = element_rect(fill = "white", color = "white"),  # Set panel background to white
          plot.background = element_rect(fill = "white", color = "white")) + 
    scale_y_continuous(limits = c(-1, 5))
  
  # Save the plots
  ggsave(paste(marker, "_regions.png", sep = ''), dpi = 300)
}
