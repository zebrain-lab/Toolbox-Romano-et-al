# Toolbox-Romano-et-al

## Computational toolbox for analysis of calcium imaging data of neuronal populations

See the preprint **A computational toolbox and step-by-step tutorial for the analysis of neuronal population dynamics in calcium imaging data** by Romano et al for tutorial and full description (doi.org/10.1101/103879).

The toolbox is demonstrated in the paper **An integrated calcium imaging processing toolbox for the analysis of neuronal population dynamics** by Romano et al (doi.org/10.1371/journal.pcbi.1005526).

Developed at Zebrafish Neuroethology lab (http://www.zebrain.biologie.ens.fr/)

## Installation
Download the latest version of the toolbox by clicking on "Clone or download". Unfortunately, GitHub does not allow downloading files larger than 50MB this way, you will obtain a compressed zip file where files >50MB will be replaced by 134-bytes "dummy" files. The test data has 3 files larger than 50MB: "OT.tif" and "OT_reg_smooth.tif" in the "Test Data" folder, and "OT_reg_smooth_RASTER.mat" in the "Test Data/Processed files" folder. To download the correct version of these large files, you should individually click on the corresponding links in the repository (e.g., navigate to "Test Data", click on "OT.tif" and then on "Download"), and replace the "dummy" 134-bytes files on your computer.


### Bug fixes
April 2017: Fixed bug in FindROIs.m that mistakenly stored inverted ROI perimeters (i.e., a "mirror image" of the correct ROIs).

### Description of variables
We now describe all the variables stored in Matlab .mat files during the utilization of the
toolbox. All the files correspond to processing and analysis of a fluorescence imaging video (e.g., myVideo.tif) with *T* imaging frames and *N* ROIs.
<br /> 
<br />
- In the **_ALL_CELLS.mat** file (e.g., myVideo_ALL_CELLS.mat), the following variables are
stored:

*avg* : average image (across frames) of the imaging file, showing the anatomy of imaged plane.

*bkg* : a logical mask containing all the ROIs found. Same size as *avg*.

*cells* : a Matlab cell array of size 1 x *N*, containing the pixel indexes of each ROI.

*cell_number* : total number of ROIs found (equal to *N*).

*cells_mean* : matrix of size *T* x *N*, containing the fluorescence time series of all the ROIs.

*cell_per* : a Matlab cell array of size *N* x 1, containing the perimeter coordinates for each ROI.

*distances* : matrix of size *N* x *N*, with the distances (in pixels) between ROIs.

*npil_mean* : matrix of size *T* x *N*, containing the fluorescence time series of the local perisomatic

neuropil that surrounds each ROIs.

*pixelLengthX* : size in micrometers of each image pixel in the X direction.

*pixelLengthY* : size in micrometers of each image pixel in the Y direction.
<br />
<br />
- In the **_ARTIFACTS.mat** file (e.g., myVideo_ARTIFACTS.mat), the following variable is stored:

*movements* : *T* x 1 binary array, with ones for imaging frames where a movement artifact was
found, and zeros otherwise.
<br />
<br />
- In the **_RASTER.mat** file (e.g., myVideo_RASTER.mat), the following variables are stored:

*raster* : a *T* x *N* matrix of the significant fluorescence transients of all accepted ROIs.

*deltaFoF* : a *T* x *N* matrix of the ∆F/F0 time series for all the accepted ROIs.

*F0* : a *T* x *N* matrix of baseline fluorescence time series matrix for all the accepted ROIs.

*deletedCells* : the indexes of the rejected ROIs.

*movements* : *T* x 1 binary array, with ones for imaging frames where a movement artifact was
detected, zero otherwise. Same variable as that in *_ARTIFACTS.mat* file.

*mu* : 1 x *N* array with the average fluorescence baseline of each accepted ROI.

*sigma* : 1 x *N* array with the estimated fluorescence baseline noise scale.

*mapOfOdds* : map of fluorescence transitions that were considered significant.

*mapOfOddsJoint* : map of fluorescence transitions that were considered significant and
biophysically realistic (rasters are determined with this map).

*xev* : x-axis of *mapOfOddsJoint*.

*yev* : y-axis of *mapOfOddsJoint*.

*densityData* : density of fluorescence transitions of all the accepted ROIs.

*densityNoise* : density of fluorescence transitions of the noise model.

*imageAvg* : the average image (across video frames) of your TIFF file.

*params* : parameters chosen by the user for the analysis of the fluorescence dynamics.
<br />
<br />
- In the **_RESPONSE_MAP.mat** file (e.g., myVideo_RESPONSE_MAP.mat), the following variables
are stored:

*traces* : Matlab cell array of size *N* x *S*, where *S* is the number of experimental events being
mapped, where the ROI ∆F/F trial traces of each trial event are stored.

*mapParams* : *S* x 1 of the parametric values of the events.

*roiHSV* : *N* x 3 matrix with original HSV color-code for the mapped ROI responses.

*roiHSVRescaled* : same as *roiHSV* after color rescaling.

*peakParameter* : *N* x 1 matrix with the event parameter that gives a ROI peak response (mapped to
hue color channel).

*peakParameterStd* : *N* x 1 matrix with the tuning width around the ROI peak response (mapped to
the inverse of the saturation color channel).

*peakResponse* : *N* x 1 matrix with the ∆F/F value of the ROI peak response (mapped to the value
color channel).
<br />
<br />
- In the **_CLUSTERS.mat** file (e.g., myVideo_CLUSTERS.mat), the following
variables are stored:

*assembliesCells* : Matlab cell array of size 1 x *M*, where *M* is the number of assemblies found,
containing the cells that participate in each assembly.

*confSynchBinary* : threshold for significance of the pooled population activity count (events of
synchronous population activity with cell count above this value are significant).

*matchIndexTimeSeries* : *M* x *T* matrix with the time series of the matching index for the
assemblies.

*matchIndexTimeSeriesSignificance* : *M* x *T* matrix with the significance (p-Value) of the time series
of the matching index for the assemblies.

*threshSignifMatchIndex* : threshold for significance used for the matching index time series of the
assemblies.

*matchIndexTimeSeriesSignificant* : *M* x *T* matrix with the complete transients of significant
assemblies' matching indexes.

*matchIndexTimeSeriesSignificantPeaks* : *M* x *T* matrix with the transient peaks of significant
assemblies' matching indexes.
<br />
<br />
- In the **_ORDER_TOPO.mat** file (e.g., myVideo_ORDER_TOPO.mat), the following variables are
stored:

*orderOfCells* : order of all ROIs, according to the projection of the assemblies centroids over the
anatomical axis drawn by the user (non-assembly cells are at the bottom).

*orderOfCellsInAssemblies* : same as *orderOfCells*, but only for assemblies' cells.

*assembliesOrdered* : Matlab cell array of size 1 x *R*, where *R* is the number of anatomical axes
selected.
<br />
<br />
- In the **_SURROGATE_CLUSTERS.mat** file (e.g., myVideo_SURROGATE_CLUSTERS.mat), the
following variables are stored:

*assembliesSurrogateRandom* : Matlab cell array of size *P* x *K*, where *K* is the number of surrogate
random assemblies per original assembly, with the list of ROIs in each random surrogate assembly.

*assembliessSurrogateTopo* : Matlab cell array of size *P* x *K*, where *K* is the number of surrogate
topographic assemblies per original assembly, with the list of ROIs in each topographic surrogate
assembly.
<br />
<br />
- In the **_ASSEMBLIES_vs_SURROGATE.mat** file (e.g., myVideo_ASSEMBLIES_vs_SURROGATE.mat), the following variables are stored:

*varAssemblies* : Matlab cell array of size 1 x *F*, where *F* is the number of features. For feature *Fi*
the variable *varAssemblies{Fi}.data* has the feature pooled according to the assemblies.

*varSurrogate* : Matlab cell array of size 1 x *F*, where *F* is the number of features. For feature *Fi* the
variables *varNMs{Fi}.Random* and *varNMs{Fi}.Topo* have the feature pooled according to the
random surrogate assemblies and the topographic surrogate assemblies, respectively.
