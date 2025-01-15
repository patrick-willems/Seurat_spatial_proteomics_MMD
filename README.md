# MMD Spatial Proteomics

Spatial proteomics image and single cell data analysis from a Moyamoya disease (MMD) biopsy.

## Spatial proteomics image analysis

### Subsetting the data

Tiff files relevant for the analysis (i.e.: excluding bleaching images, autofluorescence images, superfluous DAPI images and the first 3 cycles (which were not aligned with cycles 4 and beyond)) were copied to a local drive using the commands in copy_relevant_files_for_analysis.sh. 

The DAPI channel of cycle 4 was kept as the nucleus channel, alongside the antibody stains of 146 markers.

### Correcting images

The tiff files were corrected so that they align with each other and have the same dimensions using the notebook correct_images.ipynb. In MACS IQ this is achieved by offsetting the images using values specified in the metadata, but for analysis in sparrow, we need images with the same size and co-registered. 

After this step, the original images were removed to save space.

### Cropping images

The tiff files were subsequently cropped to a usable ROI that is in focus and does not contain fluorescent reflections or too much folded tissue. This was done using crop_images.ipynb.

After this step, the uncropped images were removed to save space.

### Region annotation

Several tissue regions (the lumen, intima, lamina elastica, tunica media and tunica adventitia) and artifacts (out-of-focus regions, tissue folds and acquisition bleaching) were annotated in Qupath and exported as a geojson (regions.geojson contains all individual regions within the image and artifacts.geojson contains all artifacts). The regions were annotated using the DAPI, Desmin and CollagenI channels. The acquisition bleaching was annotated based on the CLA channel.

### Segmentation

Segmentation was done for immune cells, RBCs separately and combined with an expanded nucleus segmentation (10 pixels expansion) using segmentation.ipynb. The gene lists used for segmentation can be found in gene_lists.py. Merging was done giving priority to immune cell segmentations (since these were the highest quality), then RBC, followed by expanded nuclei. Segmentation masks were exported as tiff.

### Registration correction of first 3 cycli

Due to problems during image acquisition, the first 3 cycli are shifted compared to the other channels in the experiment and for some reason the registration algorithm of MACSima was not able to correct this. Using manual landmarks in BigWarp, the DAPI channel of cycle 1 was registered with cycle 4. All protein channels were transformed based on the registered DAPI cycle 1 image.

The cropped registered channels of the first 3 cycli were added to the cropped_images folders where the other cropped channels are stored.

### Downstream analysis

The downsteam_analysis.ipynb was created to create an sdata object to combine all channels, segmentation masks, region annotations and artifacts annotations. Intensities were allocated to cells and corrected for cell size. Morphological cell features were calculated (e.g.: area, centroid coordinates, ...). For each cell, it was determined in which regions it occurs, whether it has a segmented nucleus or not, what the distance to the lumen is, which original segmentation group the cell belongs to (i.e. immune, RBC or expanded nucleus) and several plots were created to show these variables spatially. The tables containing the information for all cells were saved as cells_table.csv.

### Ilastik object classifiers

Object classifiers were trained in ilastik for the following immune channels: CD2, CD8a, CD15, CD43, CD45, CD45RA, CD45RB, CD57, CD66b, CD162, HLADR and Syk. CD3 was not included because the signal was too unclear.

- Object Classification (Inputs: Raw Data, Segmentation)
- Input Data: Raw Data is the cropped intensity image of the channel, Segmentation Image is the merged expanded nucleus mask.
- Object Feature Selection: Only the intensity features were used for the object classifier (every feature under Intensity Distribution)
- Object Classification: cells were classified in either positive or negative classes until the predictions looked good.
- Object Information Export: Results were exported as csv and added to cells_table.csv

## Single cell data analysis

Average intensities per cell for all channels were used to generate a 10X-like input for Seurat single cell analysis containing 138 features (antibody channels) for 3,936 cells.29 A Seurat object was generated using the Read10X function and cell metadata from the spatial proteomics image analysis was added, including the manually annotated tissue region (e.g. intima or tunica media), the presence of technical artefacts and others (Supplemental Table 1). High-quality cells were filtered by removing cells presenting outlier intensity values,acquisition bleaching or residing in out-of-focus or tissue folding areas. Further analysis for the 2,362 high-quality cells proceeded with Seurat (version 4.3.0.1) following a typical CITE-Seq analysis. For data normalization, we performed centered log ratio (CLR) transformation within cells (NormalizeData function, method ‘CLR’, margin 2). Afterwards, we detected highly variable genes (FindVariableGenes function), scaled data (ScaleData function) and performed principal component analysis (runPCA function) using 12 principal components. Subsequently, unsupervised clustering was performed using the FindNeighbours and FindClusters functions (clustering resolution of 0.2) and visualized by uniform manifold approximation and projection (UMAP) using the RunUMAP function. This reveals five clusters corresponding to distinct cell tissues (intima, adventitia, media, endothelium and lumen). Protein markers for each cluster were identified with the FindMarkers function using the ‘roc’ test (upregulated markers only), where we retained the top five ranked proteins with the highest predictive power for each cluster. This procedure was repeated for the subset of 223 immune cells identified by the whole-cell segmentation based on a combination of immune markers described above.

Script: 

## Differential protein statistics LCM LC-MS/MS

The outputted peptide precursor quantification matrix was used for differential protein analysis using MSqRob2 (version 1.6.1). We used the recommended pipeline, including log transformation of precursor intensities, requiring precursors to be quantified at least twice, and sample median normalization. Afterwards, precursor intensities were summarized at the protein level, using the aggregateFeatures function. The multidimensional scaling (MDS) plot of protein-level expression data shows VSMC and occlusion replicate samples to cluster coherently. Next, a MSqRob model was fitted using robust linear regression and the condition (VSMC or occlusion) was specified as covariate and contrast of interest. Nine differentially regulated proteins (p-value ≤ 0.01 and fold change ≥ 2) were identified.

Raw proteomics data: PRIDE PXD057183
Script: MSqRob.R

script: MACSima_analysis.R
