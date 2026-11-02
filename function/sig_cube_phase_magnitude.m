function [res_mag,res_pha]=sig_cube_phase_magnitude(w,d,ssp_map,ssm_map,Sidx_m,Sidx_p,Sv_m,Sv_p)
% sig_cube_phase_magnitude - form individual biomarker region with
% significant group difference in activation magnitude and phase polarity transitions
% Inputs:
%   - w        : cube size [w ¡Á w ¡Á w] (in voxels)
%   - d        : sampling step (in voxels)
%   - ssp_map  : source phase map within sub-region [X x Y x Z]
%            X, Y, Z represent three spatial dimensions.
%   - ssm_map  : source magnitude map within sub-region [X x Y x Z]
%   - Sidx_m   : indices of significant magnitude cubes from statistical tests
%   - Sidx_p   : indices of significant phase cubes
%   - Sv_m     : corresponding t-values for magnitude cubes
%   - Sv_p     : corresponding t-values for phase cubes
%
% Outputs:
%   - res_mag  : biomarker region of magnitude cubes with significant group difference, weighted by sign(Sv_m)
%   - res_pha  : biomarker region of magnitude cubes with significant group difference, weighted by sign(Sv_p)
%
%% Compute mean activation magnitude within magnitude voxel cubes
% calculate the number of voxel cubes extractable in each spatial dimension.
num_cubes_x = floor((size(ssp_map, 1) - w) / d) + 1;
num_cubes_y = floor((size(ssp_map, 2) - w) / d) + 1;
num_cubes_z = floor((size(ssp_map, 3) - w) / d) + 1;

% calculate mean activation magnitude
mag_cubes = zeros(num_cubes_x, num_cubes_y, num_cubes_z, w, w, w);
for nxt=1:num_cubes_x
    for nyt=1:num_cubes_y
       for nzt=1:num_cubes_z 
        mag_cubes(nxt,nyt,nzt,:,:,:)=ssm_map...
           (1+d*(nxt-1):1+w+d*(nxt-1)-1,1+d*(nyt-1):1+w+d*(nyt-1)-1,1+d*(nzt-1):1+w+d*(nzt-1)-1);
       end;end;end
mag_cubes_r=reshape(mag_cubes,[size(mag_cubes,1)*size(mag_cubes,2)*size(mag_cubes,3),w,w,w]);

for voxel_idx=1:size(mag_cubes_r,1)
   ssm=reshape(squeeze(mag_cubes_r(voxel_idx,:,:,:)),[w * w * w,1]);
   a_m(voxel_idx)=sum(ssm(:)) / ( w * w * w);
end
a_m_r=reshape(a_m,[size(mag_cubes,1),size(mag_cubes,2),size(mag_cubes,3)]);

% rank all cubes by descending magnitude
[~, sorted_ind] = sort(reshape(squeeze(a_m_r),[],1),'descend');
[order_x, order_y, order_z] = ind2sub(size(squeeze(a_m_r)), sorted_ind);


%% Reconstruct spatial map with signed t-values of mangitude
% ensure valid indexing
if max(Sidx_m)>length(order_x)
sorted_ind_x=order_x(Sidx_m(1:(min(find(Sidx_m>length(order_x)))-1)));
sorted_ind_y=order_y(Sidx_m(1:(min(find(Sidx_m>length(order_y)))-1)));
sorted_ind_z=order_z(Sidx_m(1:(min(find(Sidx_m>length(order_z)))-1)));
else 
sorted_ind_x=order_x(Sidx_m);sorted_ind_y=order_y(Sidx_m);sorted_ind_z=order_z(Sidx_m); 
end

res_mag = zeros(size(ssm_map));
for i = length(sorted_ind_x):-1:1
    cube = ssm_map(1+d*(sorted_ind_x(i)-1):1+w+d*(sorted_ind_x(i)-1)-1, ...
                   1+d*(sorted_ind_y(i)-1):1+w+d*(sorted_ind_y(i)-1)-1, ...
                   1+d*(sorted_ind_z(i)-1):1+w+d*(sorted_ind_z(i)-1)-1);
    
    t_value = sign(Sv_m(i));
    weighted_cube = cube * t_value;
    
    res_mag(1+d*(sorted_ind_x(i)-1):1+w+d*(sorted_ind_x(i)-1)-1, ...
            1+d*(sorted_ind_y(i)-1):1+w+d*(sorted_ind_y(i)-1)-1, ...
            1+d*(sorted_ind_z(i)-1):1+w+d*(sorted_ind_z(i)-1)-1) = weighted_cube;
end
clear sorted_ind_x sorted_ind_y sorted_ind_z

%% Reconstruct spatial map with signed t-values of polarity
% ensure valid indexing
if max(Sidx_p)>length(order_x)
    sorted_ind_x_ntxyz=order_x(Sidx_p(1:(min(find(Sidx_p>length(order_x)))-1)));
    sorted_ind_y_ntxyz=order_y(Sidx_p(1:(min(find(Sidx_p>length(order_y)))-1)));
    sorted_ind_z_ntxyz=order_z(Sidx_p(1:(min(find(Sidx_p>length(order_z)))-1)));
else 
    sorted_ind_x_ntxyz=order_x(Sidx_p);
    sorted_ind_y_ntxyz=order_y(Sidx_p);
    sorted_ind_z_ntxyz=order_z(Sidx_p);    
end

res_pha = zeros(size(ssp_map));
for i = length(sorted_ind_x_ntxyz):-1:1
    cube = ssm_map(1+d*(sorted_ind_x_ntxyz(i)-1):1+w+d*(sorted_ind_x_ntxyz(i)-1)-1, ...
                   1+d*(sorted_ind_y_ntxyz(i)-1):1+w+d*(sorted_ind_y_ntxyz(i)-1)-1, ...
                   1+d*(sorted_ind_z_ntxyz(i)-1):1+w+d*(sorted_ind_z_ntxyz(i)-1)-1);
    t_value = sign(Sv_p(i)); 
    weighted_cube = cube * t_value; 
    res_pha(1+d*(sorted_ind_x_ntxyz(i)-1):1+w+d*(sorted_ind_x_ntxyz(i)-1)-1, ...
            1+d*(sorted_ind_y_ntxyz(i)-1):1+w+d*(sorted_ind_y_ntxyz(i)-1)-1, ...
            1+d*(sorted_ind_z_ntxyz(i)-1):1+w+d*(sorted_ind_z_ntxyz(i)-1)-1) = weighted_cube;
end
