function [r_sm] = subregion_division(sm, com_ii, r_format)
% subregion_division - Divide a 3D spatial map into sub-regions.
%
% Inputs:
%   - sm: denoised 3D spatial map [X x Y x Z] 
%            X, Y, Z represent three spatial dimensions.
%   - com_ii: component ID:
%        1 for DMN (Default Mode Network);
%        2 for AUD (Auditory Network);
%        3 for ECN (Executive Control Network).
%   - r_format:
%        1 for returning the full 3D map for each sub-region [X x Y x Z],
%           and the areas outside the sub-region are set to zero;
%        2 for returning only the 3D sub-region.
% Output:
%   - r_sm: spatial map within sub-regions.
%%
% Define the region-specific coordinates (in voxel indices):
    % DMN
    x_acc = 16:36; y_acc = 43:63; z_acc = 13:30;   % ACC (Anterior Cingulate Cortex)    
    x_pcc = 17:37; y_pcc = 8:35; z_pcc = 21:34;    % PCC (Posterior Cingulate Cortex)
    x_ipll = 1:17; y_ipll = 1:30; z_ipll = 22:34;  % L-IPL (Left Inferior Parietal Lobule)
    x_iplr = 37:53; y_iplr = 1:30; z_iplr = 22:34; % R-IPL (Right Inferior Parietal Lobule)

    % AUD
    x_al = 1:23; y_al = 19:53; z_al = 14:26;       % L-AUD (Left Auditory)
    x_ar = 31:53; y_ar = 19:53; z_ar = 14:26;      % R-AUD (Right Auditory)

    % ECN
    x_le = 1:22; y_le = 30:63; z_le = 15:31;       % L-ECN (Left Executive Control Network)
    x_me = 22:31; y_me = 30:63; z_me = 15:31;      % M_ECN (Middle Executive Control Network)
    x_re = 31:53; y_re = 30:63; z_re = 15:31;      % R_ECN (Right Executive Control Network)

 
% Check the dimensions of the input 3D spatial map
    [sm_X, sm_Y, sm_Z] = size(sm);
    if sm_X ~= 53 || sm_Y ~= 63 || sm_Z ~= 46
        error('The input 3D spatial map dimensions should be [53 x 63 x 46].');
    end


    if r_format == 1

            if com_ii == 1
                r_sm = zeros(53, 63, 46, 4);
                r_sm(x_acc, y_acc, z_acc, 1) = sm(x_acc, y_acc, z_acc);
                r_sm(x_pcc, y_pcc, z_pcc, 2) = sm(x_pcc, y_pcc, z_pcc);
                r_sm(x_ipll, y_ipll, z_ipll, 3) = sm(x_ipll, y_ipll, z_ipll);
                r_sm(x_iplr, y_iplr, z_iplr, 4) = sm(x_iplr, y_iplr, z_iplr);
            elseif com_ii == 2
                r_sm = zeros(53, 63, 46, 2);
                r_sm(x_al, y_al, z_al, 1) = sm(x_al, y_al, z_al);
                r_sm(x_ar, y_ar, z_ar, 2) = sm(x_ar, y_ar, z_ar);
            elseif com_ii == 3
                r_sm = zeros(53, 63, 46, 3);
                r_sm(x_le, y_le, z_le, 1) = sm(x_le, y_le, z_le);
                r_sm(x_me, y_me, z_me, 2) = sm(x_me, y_me, z_me);
                r_sm(x_re, y_re, z_re, 3) = sm(x_re, y_re, z_re);
            else
                error('Invalid component index com_ii');
            end       

    else

            if com_ii == 1
                r_sm{1} = sm(x_acc, y_acc, z_acc);
                r_sm{2} = sm(x_pcc, y_pcc, z_pcc);
                r_sm{3} = sm(x_ipll, y_ipll, z_ipll);
                r_sm{4} = sm(x_iplr, y_iplr, z_iplr);
            elseif com_ii == 2
                r_sm{1} = sm(x_al, y_al, z_al);
                r_sm{2} = sm(x_ar, y_ar, z_ar);
            elseif com_ii == 3
                r_sm{1} = sm(x_le, y_le, z_le);
                r_sm{2} = sm(x_me, y_me, z_me);
                r_sm{3} = sm(x_re, y_re, z_re);
            else
               error('Invalid component index com_ii');
            end
    end
end
