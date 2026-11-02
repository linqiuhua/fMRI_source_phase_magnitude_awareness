Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase" submitted to Nature Computational Science
MATLAB Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase"
submitted to Nature Computational Science
Matlab Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase" submitted to Nature Computational Science

data  folder:
      1. inbrain_ind.mat: indices of inbrain voxels
      2. Spatial_ref.mat: spatial references (Smith et al., 2009) for DMN, AUD and ECN
      3. cvDMN_COBRE.mat: complex-valued spatial maps for DMN component from COBRE dataset (42 SZs and 40 HCs)
      4. cvDMN_UCLA.mat: complex-valued spatial maps for DMN component from UCLA dataset (50 SZs and 50 HCs)

Main code: 
       main.m : main function.

       main.m : individual-level analysis.
       GroupLevelAnalysis.m : group-level analysis.

function folder:
       1. cal_cube_phase_magnitude.m : compute the mean number of phase polarity transitions and activation magnitudes within voxel cubes.
       2. complex_fmri_scannerPhase.m : construct complex-valued fMRI data using magnitude and scanner phase data.
@@ -31,7 +32,12 @@ function folder:
       18. SM_display_multiSlice.m : multiple slices visualization of a spatial map.
       19. subregion_division.m : divide a spatial map into sub-regions.
       20. toolbox : functions (from GIFT toolbox) for NIfTI data and ICA analysis.

       21. Boxplot.m : generate boxplots of the piecewise linear model parameters and fitting errors for healthy controls and schizophrenia patients, marking significant group differences.
       22. FDR_group_differences.m : perform group-level difference with BH-FDR correction for activation magnitude and phase polarity transitions.
       23. sig_cube_phase_magnitude.m : form individual biomarker region with significant group-level difference in activation magnitude and phase polarity transitions.
       24. SMdisplay_overlap_fdr.m : 3D visualization of the overlapping spatial maps between significant difference of activation magnitude and phase polarity transitions after BH-FDR correction.
       25. SVM_classification.m : perform SVM classification using using activation magnitude, phase polarity transitions and their fusion. 

SMshow folder: 
       code and utilities for visualizing 3D brain maps

@@ -41,3 +47,4 @@ Environment:
Data:
       Complex-valued fMRI dataset: https://doi.org/10.5281/zenodo.20619435
       Magnitude-onlye fMRI dataset: https://openneuro.org/datasets/ds000030
