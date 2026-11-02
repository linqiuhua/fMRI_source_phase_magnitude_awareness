function [n_t, a_m] = cal_cube_phase_magnitude(w, d, ssp_map, ssm_map)
% cal_cube_phase_magnitude - Compute the mean phase polarity transitions and mean activation magnitude 
%                            within local voxel cubes, based on phase and magnitude maps respectively.

% Inputs:
%   - w: cube size [w x w x w] (in voxels).
%   - d: sampling step (in voxels).
%   - ssp_map: denoised phase map [X x Y x Z].  X, Y, Z represent three spatial dimensions.
%   - ssm_map: denoised magnitude map [X x Y x Z].
%
% Outputs:
%   - n_t:  mean number of phase polarity transitions along x, y, z directions 
%             within each phase voxel cube [3xN]
%             N = number of phase voxel cubes
%   - a_m: mean activation magnitude within each magnitude voxel cube [1xN]
%             N = number of magnitude voxel cubes

%% Phase Analysis: Calculate Phase Polarity Transitions
%-- Step 1: phase cube sampling --%
    % Calculate the number of cubes along each axis
    num_cubes_x = floor((size(ssp_map, 1) - w) / d) + 1;
    num_cubes_y = floor((size(ssp_map, 2) - w) / d) + 1;
    num_cubes_z = floor((size(ssp_map, 3) - w) / d) + 1;

    phase_cubes = zeros(num_cubes_x, num_cubes_y, num_cubes_z, w, w, w);
    for idx_x = 1:num_cubes_x
        for idx_y = 1:num_cubes_y
            for idx_z = 1:num_cubes_z
                phase_cubes(idx_x, idx_y, idx_z, :, :, :) = ssp_map(1 + d * (idx_x - 1):w + d * (idx_x - 1), ...
                                                               1 + d * (idx_y - 1):w + d * (idx_y - 1), ...
                                                               1 + d * (idx_z - 1):w + d * (idx_z - 1));
            end
        end
    end

    phase_cubes = reshape(phase_cubes, [size(phase_cubes, 1) * size(phase_cubes, 2) * size(phase_cubes, 3), w, w, w]);

%-- Step 2: calculate mean phase polarity transitions --%
    n_t = zeros(3, size(phase_cubes, 1));     
    
    for voxel_idx = 1:size(phase_cubes, 1)
    % x direction
            nt_x = zeros(w, w);    
            for i = 1:w
                for j = 1:w
                    phase_cube_x = phase_cubes(voxel_idx, :, i, j);
                    % Count polarity transitions (change of sign)
                    nt_x(i, j) = sum(phase_cube_x(1:w-1) .* phase_cube_x(2:w) < 0);
                end
            end

           % Calculate the mean of non-zero transitions
            nonzero_transitions_x = nt_x(nt_x ~= 0);  

            if isempty(nonzero_transitions_x)
                n_t(1, voxel_idx) = 0;   % If no non-zero transitions, set to 0
            else
                n_t(1, voxel_idx) = mean(nonzero_transitions_x);
            end    
    
    % y direction
            nt_y = zeros(w, w);    
            for i = 1:w
                for j = 1:w
                    phase_cube_y = phase_cubes(voxel_idx, i, :, j); 
                    nt_y(i, j) = sum(phase_cube_y(1:w-1) .* phase_cube_y(2:w) < 0);
                end
            end
            nonzero_transitions_y = nt_y(nt_y ~= 0);  
            if isempty(nonzero_transitions_y)
                n_t(2, voxel_idx) = 0;
            else
                n_t(2, voxel_idx) = mean(nonzero_transitions_y);
            end
    
       % z direction
            nt_z = zeros(w, w);    
            for i = 1:w
                for j = 1:w
                    phase_cube_z = phase_cubes(voxel_idx, i, j, :); 
                    nt_z(i, j) = sum(phase_cube_z(1:w-1) .* phase_cube_z(2:w) < 0);
                end
            end

            nonzero_transitions_z = nt_z(nt_z ~= 0);  
            if isempty(nonzero_transitions_z)
                n_t(3, voxel_idx) = 0;
            else
                n_t(3, voxel_idx) = mean(nonzero_transitions_z);
            end
end   
    
%% Magnitude Analysis: Calculate Activation Magnitude
%-- Step 1: magnitude cube sampling --%
    mag_cubes = zeros(num_cubes_x, num_cubes_y, num_cubes_z, w, w, w);
    for idx_x = 1:num_cubes_x
        for idx_y = 1:num_cubes_y
            for idx_z = 1:num_cubes_z
                mag_cubes(idx_x, idx_y, idx_z, :, :, :) = ssm_map(1 + d * (idx_x - 1):w + d * (idx_x - 1), ...
                                                                   1 + d * (idx_y - 1):w + d * (idx_y - 1), ...
                                                                   1 + d * (idx_z - 1):w + d * (idx_z - 1));
            end
        end
    end
    
    mag_cubes = reshape(mag_cubes, [size(mag_cubes, 1) * size(mag_cubes, 2) * size(mag_cubes, 3), w, w, w]);

%-- Step 2: calculate mean activation magnitude --%
    a_m = zeros(1, size(mag_cubes, 1));

    for voxel_idx = 1:size(mag_cubes, 1)
        ssm = reshape(squeeze(mag_cubes(voxel_idx, :, :, :)), [w * w * w, 1]);
        a_m(voxel_idx) = sum(ssm(:)) / (w * w * w); 
    end

end
