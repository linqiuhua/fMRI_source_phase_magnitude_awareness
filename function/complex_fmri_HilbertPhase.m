function data_complex_inbrain = complex_fmri_HilbertPhase(mag_file, phase_data,inbrain_ind)
% complex_fmri_construct - Construct complex-valued fMRI data from magnitude NIfTI and Hilbert-transform generated phase.
%
%   Inputs:
%    - mag_file: file path to the magnitude NIfTI image.
%    - phase_data: Phase data derived from Hilbert transform.[V x 1].
%    - inbrain_ind: inbrain voxel indices [V x 1]
%
%   Output:
%    - data_complex_inbrain: inbrain complex-valued fMRI data [T x V] 
%                            T = time points, V = number of inbrain voxels.


%%
% Load magnitude data
mag_nii = load_nii(mag_file);
mag_data = double(mag_nii.img(:,:,:,1:146));  % Exclude unstable timepoints
[nx, ny, nz, T] = size(mag_data);
mag_2d = reshape(mag_data, [nx * ny * nz, T]).';
mag_inbrain = mag_2d(:, inbrain_ind);

i = sqrt(-1);
data_complex_inbrain = mag_inbrain.*exp(i*phase_data);
end
