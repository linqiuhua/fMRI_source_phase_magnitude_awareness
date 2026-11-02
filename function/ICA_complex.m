function [SS_all, At_all] = ICA_complex(data, ModelOrder, Num_run)
%ICA_complex - Perform complex ICA on fMRI data across multiple model orders and runs.
%
%   Inputs:
%    - data: inbrain complex-valued fMRI data [T x V],
%             T = time points, V = number of in-brain voxels.
%    - ModelOrder: vector of component numbers, e.g., [20, 40, 60].
%    - Num_run: number of ICA separation per model order.
%
%   Outputs:
%    - SS_all: cell array of spatial components [Num_run x Num_order],
%              each cell contains a [MO x V] matrix.
%    - At_all: cell array of time courses [Num_run x Num_order],
%              each cell contains a [T x MO] matrix.
%% Initialize
Num_order = length(ModelOrder);
SS_all = cell(Num_run, Num_order);
At_all = cell(Num_run, Num_order);

%% ICA loop over runs and model orders
parfor r = 1 : Num_run
   % Remove mean per time point
    data_pro = icatb_preproc_data(data.', 'remove mean per timepoint');
    
    for m = 1: Num_order
        % Principal component analysis 
        [data_pca, ~, ~, ~, ~] = ...
              icatb_calculate_pca(data_pro, ModelOrder(m), 'remove_mean', 0); 
        % Complex-valued ICA decomposition
        [W, ~, Se] = complex_ICA_EBM(data_pca.'); 
        % Time courses estimation
        At_all {r,m} = data * pinv(Se); 
        % Spatial component reconstruction
        SS_all {r,m} = W * data_pca.'; 
    end
end
