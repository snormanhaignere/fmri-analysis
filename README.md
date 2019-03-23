# fmri-analysis
General purpose code for fMRI analysis

This code is provided mainly as a reference for our publication entitled: 

Divergence in the Functional Organization of Human and Macaque Auditory Cortex Revealed by fMRI Responses to Harmonic Tones

Some of the key functions are described below. The scripts are not intended to be stand alone (and thus will not run unless files are formatted appropriately), but nonetheless provide a useful reference for the analyses described in the paper. The repository fmri-analysis-v2 provides a more up-to-date, user-friendly pipeline.

All of the analyses were performed on 

fla_matlab.m: First level GLM analyses and corresponding permutation tests.

tla_permutation_voxelwise_stats.m: Second or third level analyses combining across runs from a subject (or across subjects, but this is not relevant for the paper noted above). This script computes voxel-wise p-values using a permutation test across stimuli.

tla_permutation_clustcorr.m: Cluster-corrected maps also based on a permutation test. Assumes voxel-wise stats have already been computed.

roi_surf_grid_v2.m: ROI analyses.

downsample_surface_timecourses.m: Downsamples the finely-sampled surface mesh computed by Freesurfer to a 2D surface grid that easier/faster to work with.