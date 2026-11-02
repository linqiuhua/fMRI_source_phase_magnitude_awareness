function data_phase_inbrain = phase_fmri_Hilbert(mag_file, inbrain_ind)
% phase_fmri_Hilbert - Generate fMRI phase data by applying the Hilbert transform to the magnitude NIfTI.
%
%   Inputs:
%    - mag_file: file path to the magnitude NIfTI image.
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
data_complex = hilbert_transform(mag_inbrain);
data_phase_inbrain=angle(data_complex);
end
