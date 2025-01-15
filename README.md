# MMD Spatial Proteomics

Seurat single cell analysis based on spatial proteomics data from a Moyamoya disease (MMD) biopsy.

## Subsetting the data

Tiff files relevant for the analysis (i.e.: excluding bleaching images, autofluorescence images, superfluous DAPI images and the first 3 cycles (which were not aligned with cycles 4 and beyond)) were copied to a local drive using the commands in copy_relevant_files_for_analysis.sh. 

The DAPI channel of cycle 4 was kept as the nucleus channel, alongside the antibody stains of 146 markers.

## Correcting images

The tiff files were corrected so that they align with each other and have the same dimensions using the notebook correct_images.ipynb. In MACS IQ this is achieved by offsetting the images using values specified in the metadata, but for analysis in sparrow, we need images with the same size and co-registered. 

After this step, the original images were removed to save space.

## Cropping images

The tiff files were subsequently cropped to a usable ROI that is in focus and does not contain fluorescent reflections or too much folded tissue. This was done using crop_images.ipynb.

After this step, the uncropped images were removed to save space.

## Region annotation

Several tissue regions (the lumen, intima, lamina elastica, tunica media and tunica adventitia) and artifacts (out-of-focus regions, tissue folds and acquisition bleaching) were annotated in Qupath and exported as a geojson (regions.geojson contains all individual regions within the image and artifacts.geojson contains all artifacts). The regions were annotated using the DAPI, Desmin and CollagenI channels. The acquisition bleaching was annotated based on the CLA channel.

## Segmentation

Segmentation was done for immune cells, RBCs separately and combined with an expanded nucleus segmentation (10 pixels expansion) using segmentation.ipynb. The gene lists used for segmentation can be found in gene_lists.py. Merging was done giving priority to immune cell segmentations (since these were the highest quality), then RBC, followed by expanded nuclei. Segmentation masks were exported as tiff.

## Registration correction of first 3 cycli

Due to problems during image acquisition, the first 3 cycli are shifted compared to the other channels in the experiment and for some reason the registration algorithm of MACSima was not able to correct this. Using manual landmarks in BigWarp, the DAPI channel of cycle 1 was registered with cycle 4. All protein channels were transformed based on the registered DAPI cycle 1 image.

The cropped registered channels of the first 3 cycli were added to the cropped_images folders where the other cropped channels are stored.

## Downstream analysis

The downsteam_analysis.ipynb was created to create an sdata object to combine all channels, segmentation masks, region annotations and artifacts annotations. Intensities were allocated to cells and corrected for cell size. Morphological cell features were calculated (e.g.: area, centroid coordinates, ...). For each cell, it was determined in which regions it occurs, whether it has a segmented nucleus or not, what the distance to the lumen is, which original segmentation group the cell belongs to (i.e. immune, RBC or expanded nucleus) and several plots were created to show these variables spatially. The tables containing the information for all cells were saved as cells_table.csv.

## Ilastik object classifiers

Object classifiers were trained in ilastik for the following immune channels: CD2, CD8a, CD15, CD43, CD45, CD45RA, CD45RB, CD57, CD66b, CD162, HLADR and Syk. CD3 was not included because the signal was too unclear.

- Object Classification (Inputs: Raw Data, Segmentation)
- Input Data: Raw Data is the cropped intensity image of the channel, Segmentation Image is the merged expanded nucleus mask.
- Object Feature Selection: Only the intensity features were used for the object classifier (every feature under Intensity Distribution)
- Object Classification: cells were classified in either positive or negative classes until the predictions looked good.
- Object Information Export: Results were exported as csv and added to cells_table.csv

