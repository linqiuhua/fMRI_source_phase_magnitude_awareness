# fMRI_source_phase_magnitude_awareness
Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase" submitted to Nature Computational Science
MATLAB Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase"
submitted to Nature Computational Science

data  folder:
      1. inbrain_ind.mat: indices of inbrain voxels
      2. Spatial_ref.mat: spatial references (Smith et al., 2009) for DMN, AUD and ECN
      3. sm_COBRE82subjects_DMN.mat: complex-valued spatial maps for DMN component from 40 healthy controls and 42 schizophrenia patients in COBRE dataset
      4. sm_UCLA100subjects_DMN.mat: complex-valued spatial maps for DMN component from 50 healthy controls and 50 schizophrenia patients in UCLA dataset

Main code: 
       main.m : individual-level analysis.
       Grouplevel_analysis.m : group-level analysis.
       
function folder:
       1. Boxplot.m : Generates boxplots of the piecewise linear model parameters and fitting errors for healthy controls and schizophrenia patients, marking significant group differences.
       2. cal_cube_phase_magnitude.m : compute the mean number of phase polarity transitions and activation magnitudes within voxel cubes.
       3. complex_fmri_scannerPhase.m : construct complex-valued fMRI data using magnitude and scanner phase data.
       4. complex_fmri_HilbertPhase.m : construct complex-valued fMRI data using magnitude and Hilbert-transform phase data.
       5. complex_ICA_EBM.m : perform complex-valued ICA using EBM algorithm to extract spatial source magnitude and phase.
       6. disp_PeakCurve.m :  display magnitude and phase polarity curves around the region peak voxel.
       7. energymax_T1.m : perform phase de-ambiguity for a spatial map estimated by ICA.
       8. FDR_group_differences.m : perform group-level comparisons with BH-FDR correction for activation magnitude and phase polarity transitions.
       9. hilbert_transform.m : compute instantaneous phase via Hilbert transform.
       10. ICA_complex.m : perform complex-valued ICA (EBM).
       11. ICA_real.m : perform real-valued ICA (Infomax). 
       12. Infomax_Options.mat : Infomax ICA parameter configurations.
       13. phase_denoising.m ：perform phase de-noising (±pi/4).
       14. phase_fmri_Hilbert.m ：generate fMRI phase data by applying the Hilbert transform to the magnitude-only data.
       15. piecewise_linear_fit.m : perform adaptive piecewise linear regression to fit activation magnitudes and phase polarity transitions.
       16. remove_outliers.m : remove statistical outliers using the interquartile range rule (IQR) before fitting.
       17. selection_component_complex.m : select complex-valued component of interest.
       18. selection_component_real.m : select magnitude-only component of interest.
       19. sig_cube_phase_magnitude.m : form individual biomarker region with significant group difference in activation magnitude and phase polarity transitions.
       20. SMdisplay.m : 3D visualization of a spatial map.
       21. SMdisplay_overlap_fdr.m 3D visualization of the overlapping spatial maps between significant difference of activation magnitude and phase polarity transitions after BH-FDR correction.
       22. SM_display_multiSlice.m : multiple slices visualization of a spatial map.
       23. subregion_division.m : divide a spatial map into sub-regions.
       24. SVM_classification.m : perform SVM classification. 
       25. toolbox : functions (from GIFT toolbox) for NIfTI data and ICA analysis.

SMshow folder: 
       code and utilities for visualizing 3D brain maps
       
Environment:
       MATLAB R2014a

Data:
       Complex-valued COBRE fMRI dataset: https://doi.org/10.5281/zenodo.20619435
       Magnitude-only UCLA fMRI dataset: https://openneuro.org/datasets/ds000030
