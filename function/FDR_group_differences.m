
function result = FDR_group_differences( ...
    pvalue_am, pvalue_ntxyz, order_am, common_ind, ...
     alpha, valid_start, valid_end, ...
    tvalue_am, tvalue_ntxyz)

%% ============================================================
% Group differences and BH-FDR correction
%
% Description:
%   Perform separate BH-FDR corrections on the p-values of
%   amplitude and phase-transition features, respectively.
%
% Inputs:
%   pvalue_am       : Original p-values of amplitude features
%   pvalue_ntxyz    : Original p-values of phase-transition features
%   order_am        : Data used to determine the total number of cubes
%   common_ind      : Common valid voxel indices
%   method          : FDR correction method ('separate')
%   alpha           : Significance level (default: 0.05)
%   valid_start     : Starting index of the valid cube range
%   valid_end       : Ending index of the valid cube range
%   tvalue_am       : t-values of amplitude features
%   tvalue_ntxyz    : t-values of phase-transition features
%   do_plot         : Whether to generate plots
%
% Output:
%   result          : Structure containing FDR-corrected p-values,
%                     significant cube indices, and summary statistics
%
% Dependency:
%   mafdr (Bioinformatics Toolbox)
%% ============================================================


%% -------------------- Parameter initialization --------------------

nCube = size(order_am, 2);

if isempty(alpha)
    alpha = 0.05;
end

if isempty(valid_start)
    valid_start = 1;
end

if isempty(valid_end)
    valid_end = nCube - max(common_ind);
end

valid_end = min(valid_end, nCube);
valid_idx = valid_start:valid_end;

% Convert p-values into row vectors
pvalue_am = pvalue_am(:)';
pvalue_ntxyz = pvalue_ntxyz(:)';


%% -------------------- Extract valid p-values --------------------

% Retain p-values within the specified valid cube range
p_am_valid = pvalue_am(valid_idx);
p_ntxyz_valid = pvalue_ntxyz(valid_idx);

% Exclude NaN values from FDR correction
valid_am_mask = ~isnan(p_am_valid);
valid_ntxyz_mask = ~isnan(p_ntxyz_valid);

valid_idx_am = valid_idx(valid_am_mask);
valid_idx_ntxyz = valid_idx(valid_ntxyz_mask);

p_am_valid = p_am_valid(valid_am_mask);
p_ntxyz_valid = p_ntxyz_valid(valid_ntxyz_mask);


%% -------------------- Separate BH-FDR correction --------------------

% Apply BH-FDR correction independently to the two feature types
p_fdr_am_valid = run_bhfdr(p_am_valid);
p_fdr_ntxyz_valid = run_bhfdr(p_ntxyz_valid);


%% -------------------- Restore corrected p-values --------------------

% Initialize full-length arrays with NaN
p_fdr_am = nan(1, nCube);
p_fdr_ntxyz = nan(1, nCube);

% Map corrected p-values back to their original cube indices
p_fdr_am(valid_idx_am) = p_fdr_am_valid;
p_fdr_ntxyz(valid_idx_ntxyz) = p_fdr_ntxyz_valid;


%% -------------------- Identify significant cubes --------------------

% Significant cubes before FDR correction
sig_am_uncorrected = find( ...
    pvalue_am < alpha & ismember(1:nCube, valid_idx));

sig_ntxyz_uncorrected = find( ...
    pvalue_ntxyz < alpha & ismember(1:nCube, valid_idx));

% Significant cubes after BH-FDR correction
sig_am_fdr = find(p_fdr_am < alpha);
sig_ntxyz_fdr = find(p_fdr_ntxyz < alpha);


%% -------------------- Store analysis results --------------------

result = struct();

result.method = 'separate';
result.alpha = alpha;
result.nCube = nCube;
result.valid_idx = valid_idx;

% FDR-corrected p-values
result.p_fdr_am = p_fdr_am;
result.p_fdr_ntxyz = p_fdr_ntxyz;

% Significant cube indices
result.sig_am_uncorrected = sig_am_uncorrected;
result.sig_ntxyz_uncorrected = sig_ntxyz_uncorrected;

result.sig_am_fdr = sig_am_fdr;
result.sig_ntxyz_fdr = sig_ntxyz_fdr;

% Number of significant cubes
result.n_sig_am_uncorrected = length(sig_am_uncorrected);
result.n_sig_ntxyz_uncorrected = length(sig_ntxyz_uncorrected);

result.n_sig_am_fdr = length(sig_am_fdr);
result.n_sig_ntxyz_fdr = length(sig_ntxyz_fdr);

% Store t-values
result.tvalue_am = tvalue_am;
result.tvalue_ntxyz = tvalue_ntxyz;

result.figure_handle = [];
result.axes_handle = [];

end

%% ============================================================
% Subfunction: Benjamini-Hochberg FDR correction
%% ============================================================

function p_fdr = run_bhfdr(pvalue)

if isempty(pvalue)
    p_fdr = [];
else
    p_fdr = mafdr(pvalue(:), 'BHFDR', true);
end

end
