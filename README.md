# fMRI_source_phase_magnitude_awareness
Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase" submitted to Nature Computational Science
MATLAB Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase"
submitted to Nature Computational Science

data  folder:
      1. inbrain_ind.mat: indices of inbrain voxels
      2. Spatial_ref.mat: spatial references (Smith et al., 2009) for DMN, AUD and ECN

Main code: 
       main.m : main function.
       
function folder:
       1. cal_cube_phase_magnitude.m : compute the mean number of phase polarity transitions and activation magnitudes within voxel cubes.
       2. complex_fmri_scannerPhase.m : construct complex-valued fMRI data using magnitude and scanner phase data.
       3. complex_fmri_HilbertPhase.m : construct complex-valued fMRI data using magnitude and Hilbert-transform phase data.
       4. complex_ICA_EBM.m : perform complex-valued ICA using EBM algrithm to extract spatial source magnitude and phase.
       5. disp_PeakCurve.m:  display magnitude and phase polarity curves around the region peak voxel.
       6. energymax_T1.m : perform phase de-ambiguity for a spatial map estimated bt ICA.
       7. hilbert_transform.m : compute instantaneous phase via Hilbert transform.
       8. ICA_complex.m : perform complex-valued ICA (EBM).
       9. ICA_real.m : perform real-valued ICA (Infomax). 
       10. Infomax_Options.mat : Infomax ICA parameter configurations.
       11. phase_denoising.m ：perform phase de-noising (±pi/4).
       12. phase_fmri_Hilbert.m ：generate fMRI phase data by applying the Hilbert transform to the magnitude-only data.
       13. piecewise_linear_fit.m : perform adaptive piecewise linear regression to fit activation magnitudes and phase polarity transitions.
       14. remove_outliers.m : remove statistical outliers using the interquartile range rule (IQR) before fitting.
       15. selection_component_complex.m : select complex-valued component of interest.
       16. selection_component_real.m : select magnitude-only component of interest.
       17. SMdisplay.m : 3D visualization of a spatial map.
       18. SM_display_multiSlice.m : multiple slices visualization of a spatial map.
       19. subregion_division.m : divide a spatial map into sub-regions.
       20. toolbox : functions (from GIFT toolbox) for NIfTI data and ICA analysis.

SMshow folder: 
       code and utilities for visualizing 3D brain maps

example_results folder:
       exampling results of running main.m
       
Environment:
       MATLAB R2014a

Data:
       Complex-valued fMRI dataset: https://doi.org/10.5281/zenodo.20619435
       Magnitude-onlye fMRI dataset: https://openneuro.org/datasets/ds000030
