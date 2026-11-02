function data_complex_inbrain = complex_fmri_scannerPhase(mag_file, phase_file, inbrain_ind)
%complex_fmri_real - Generate complex-valued fMRI data using magnitude and phase NIfTI inputs.
%
%   Inputs:
%    - mag_file: file path to the magnitude NIfTI image.
%    - phase_file: file path to the phase NIfTI image. Required for 'RealComplex' method.
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

% Check if phase file is provided
if isempty(phase_file) || exist(phase_file, 'file') ~= 2
  error('Phase file must be provided and exist for ''RealComplex'' method.');
end
% Load phase data
phase_nii = load_nii(phase_file);
phase_data = double(phase_nii.img(:,:,:,1:146));
phase_2d = reshape(phase_data, [nx * ny * nz, T]).';

% Construct real complex-valued data using magnitude and phase
data_complex = mag_2d .* exp(1i * phase_2d);
data_complex_inbrain = data_complex(:, inbrain_ind);  % T ¡Á N

end
