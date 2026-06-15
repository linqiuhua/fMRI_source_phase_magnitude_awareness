%**************************************************************************
% Description: Code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase"
% submitted to Nature Computational Science
%
% Functionality:
%   - Construct complex-valued fMRI data for a single subject, phase data from either scanner or Hilbert transform
%   - Perform complex-valued ICA (EBM) to extract spatial source phase and magnitude maps
%   - Perform real-valued ICA (Infomax) to extract spatial magnitude maps
%   - Visualize spatial source phase maps and polarity curves passing through the peak voxel
%   - Perform piecewise linear modeling based on phase polarity transitions and activation magnitude
%
% Input:
%   - Preprocessed fMRI magnitude data (4D NIFTI, [X x Y x Z x T]) 
%              Optional scanner phase data (if available)
%   - inbrain_ind.mat: indices of in-brain voxels [1 x V]
%              V = number of in-brain voxels
%   - Spatial_ref.mat: spaital references (Smith et al., 2009) for three brain functional networks of interest
%              1 - DMN (Default Mode Network)
%              2 - AUD (Auditory Network)
%              3 - ECN (Executive Control Network)
%      
% Output:
%    - Spatial source phase and magnitude maps (complex-valued data) or magnitude maps (magnitude-only data)
%    - Phase polarity curves around the peak activation voxel (complex-valued data)
%    - Scatter plot and piecewise linear regression between phase polarity transitions and activation magnitude
%    - Params: structure containing model fiting parameters:
%          .b0 : intercept
%          .b1 : slope of first segment
%          .b2 : slope of second segment
%          .a_b: breakpoint of activation magnitude
%          .n_b: breakpoint of phase polarity transitions
%   - error: fitting residuals of the piecewise model
%**************************************************************************

clear; close all
%% ----------------------------Add path for self-defined functions----------------------- %%  
hSelect = helpdlg('Please select the "function" subfolder.', 'Select Function Folder');
uiwait(hSelect);
folderPath = uigetdir(pwd);

if folderPath == 0
    disp('No folder selected. Exiting.');
    return;
end
addpath(folderPath);
folderContents = dir(folderPath);
for i = 1:length(folderContents)
    if folderContents(i).isdir && ~strcmp(folderContents(i).name, '.') && ~strcmp(folderContents(i).name, '..')
        subfolderPath = fullfile(folderPath, folderContents(i).name);
        addpath(genpath(subfolderPath));
    end
end
savepath;

%% ---------- Selection of fMRI data (either magnitude-only or complex-valued) ---------- %%
fprintf('\n==================== fMRI data analysis ====================\n');
fprintf('Please select the type of fMRI analysis:\n');
fprintf('  1 - Real-valued fMRI analysis\n');
fprintf('  2 - Complex-valued fMRI analysis\n');
fprintf('==============================================================\n');

analysisMode = input('Enter 1 or 2: ');
while ~ismember(analysisMode, [1, 2])
    analysisMode = input('Invalid selection. Please enter 1 (Magnitude-only fMRI data) or 2 (Complex-valued fMRI data): ');
end

%% --------- Construct complex-valued fMRI data or load magnitude-only fMRI data ---------- %%
if analysisMode == 2
    dataType = input([ ...
        'Please select the mode of complex-valued input data:\n' ...
        '  1 - Magnitude + Scanner phase\n' ...
        '  2 - Magnitude only + Hilbert-transform generated phase \n' ...
        'Enter 1 or 2: ']);
    while ~ismember(dataType, [1, 2])
        dataType = input('Invalid selection. Please enter 1 or 2: ');
    end
else
    fprintf('Magnitude-only data is required.\n');
    dataType = 2;
end

% Load in-brain voxel indices
load('data/inbrain_ind.mat');

% Define default folder for data selection
defaultDataFolder = fullfile(pwd, 'data');
if ~exist(defaultDataFolder, 'dir')
    warning('Default data folder does not exist. Using current directory instead.');
    defaultDataFolder = pwd;
end

if analysisMode == 2 && dataType == 1
    fprintf('Input mode: Magnitude + Phase. \n Constructing complex-valued fMRI data...\n');
elseif analysisMode == 2 && dataType == 2
    fprintf('Input mode: Magnitude only. \n Estimating phase via Hilbert transform and constructing complex-valued data...\n');
else
    fprintf('Input mode: Magnitude only (Run real-valued ICA Infomax).\n');
end

switch dataType
    case 1  % Magnitude + Phase input
        hMag = helpdlg('Please select the Magnitude NIfTI data.', 'Select Magnitude Data');
        uiwait(hMag);
        [magFile, magPath] = uigetfile('*.nii', 'Select Magnitude NIfTI Data', defaultDataFolder);
        if isequal(magFile, 0); error('No magnitude data selected.'); end
        magnitudeFilePath = fullfile(magPath, magFile);

        hPhase = helpdlg('Please select the Phase NIfTI data.', 'Select Phase Data');
        uiwait(hPhase);
        [phaseFile, phasePath] = uigetfile('*.nii', 'Select Phase NIfTI Data', defaultDataFolder);
        if isequal(phaseFile, 0); error('No phase data selected.'); end
        phaseFilePath = fullfile(phasePath, phaseFile);
        
        % Generate complex-valued fMRI data using magnitude and phase
        data_for_ica = complex_fmri_scannerPhase(magnitudeFilePath, phaseFilePath, inbrain_ind);

    case 2
        hMag = helpdlg('Please select the Magnitude NIfTI data.', 'Select Magnitude Data');
        uiwait(hMag);
        [magFile, magPath] = uigetfile('*.nii', 'Select Magnitude NIfTI Data', defaultDataFolder);
        if isequal(magFile, 0); error('No magnitude data selected.'); end
        magnitudeFilePath = fullfile(magPath, magFile);
        
        if analysisMode == 2 % Complex-valued data analysis (magnitude-only + Hilbert transform phase)
            H_phase = phase_fmri_Hilbert(magnitudeFilePath, inbrain_ind);
            data_for_ica = complex_fmri_HilbertPhase(magnitudeFilePath, H_phase, inbrain_ind);
        
        else % Magnitude-only data
            nii = load_nii(magnitudeFilePath);
            magnitude_4D = nii.img;
            [X, Y, Z, T] = size(magnitude_4D);
            magnitude_2D = reshape(magnitude_4D, X*Y*Z, T);
            data_for_ica = magnitude_2D(inbrain_ind, :); 
        end
end

%% ------ Extraction of source magnitude and phase maps ------- %%
%--- Step1: Perform ICA (real or complex) ---%
% Set default ICA parameters
defaultModelOrders = [40, 50, 60]; % Model orders
defaultNumRuns = 4;                % Number of ICA runs

modelOrdersInput = input(['Enter ICA model order(s) as a vector (e.g. [40,50,60])\n' ...
                         '(press Enter for default: [' num2str(defaultModelOrders) ']): ']);
if isempty(modelOrdersInput)
    modelOrders = defaultModelOrders;
else
    modelOrders = modelOrdersInput;
end

numRunsInput = input(['Enter number of ICA runs\n' ...
                     '(press Enter for default: ' num2str(defaultNumRuns) '): ']);
if isempty(numRunsInput)
    numRuns = defaultNumRuns;
else
    numRuns = numRunsInput;
end

% Run ICA decomposition
if analysisMode == 1
    fprintf('Running real-valued ICA (model order(s): %s, runs: %d)...\n', ...
            mat2str(modelOrders), numRuns);
   
    [Sm, Tc] = ICA_real(data_for_ica, modelOrders, numRuns);
    icaTypeStr = 'real';
else
    fprintf('Running complex-valued ICA (model order(s): %s, runs: %d)...\n', ...
            mat2str(modelOrders), numRuns);
    [Sm, Tc] = ICA_complex(data_for_ica, modelOrders, numRuns);
    icaTypeStr = 'complex';
end

% Save results
if analysisMode == 1
    dataTypeStr = 'RealICA';
else
    if dataType == 1
        dataTypeStr = 'ComplexICA-ScannerPhase';
    else
        dataTypeStr = 'ComplexICA-HilbertPhase';
    end
end

resultsFolder = fullfile(pwd, 'example_results');
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

[~, personName, ~] = fileparts(magFile);
pattern = sprintf('%s_ID(\\d+)_', personName);
existingFiles = dir(fullfile(resultsFolder, sprintf('%s_ID*_*.mat', personName)));
maxID = 0;
for k = 1:numel(existingFiles)
    tok = regexp(existingFiles(k).name, pattern, 'tokens', 'once');
    if ~isempty(tok)
        idNum = str2double(tok{1});
        if ~isnan(idNum) && idNum > maxID
            maxID = idNum;
        end
    end
end
infoID = maxID + 1;

modelOrdersStr = strjoin(arrayfun(@num2str, modelOrders, 'UniformOutput', false), '_');
runStr = sprintf('%druns', numRuns);
dtStr = datestr(now, 'yyyymmdd_HHMMSS');

fileName = sprintf('%s_ID%03d_%s_%s_mo%s_%s_%s.mat', ...
                   personName, infoID, dataTypeStr, icaTypeStr, modelOrdersStr, runStr, dtStr);
fprintf('Saving ICA results to: %s\n', fullfile(resultsFolder, fileName));

save(fullfile(resultsFolder, fileName), ...
     'Sm', 'Tc', 'modelOrders', 'numRuns', 'magnitudeFilePath', 'analysisMode');

%--- Step2: Select component of interest ---%
% Load spatial reference maps (e.g., DMN, AUD, ECN)
load('data/Spatial_ref.mat');

% Select component
fprintf('Extracting components of interest...\n');

if analysisMode == 2
    for num_com = 1 : 3
        Sref = Sr(num_com,:);
        Sref(Sref < 2.5) = 0;
        [SM_select(:,:,num_com,:), SM_best(:,num_com,:)] = selection_component_complex(Sm, Tc, Sref);
    end
else
    for num_com = 1 : 3
        Sref = Sr(num_com,:);
        Sref(Sref < 2.5) = 0;       
        [SM_select(:,:,num_com,:), SM_best(:,num_com,:)] = selection_component_real(Sm, Tc, Sref);       
    end
end

%% ----------------- Spatial visualization + Fitting modeling ----------------- %%
continueFlag = true;

while continueFlag

    %-- Step 1: select component and sub-region ---%
    com_ID = input(['Select a component to visualize:\n' ...
                       '  1 - DMN (Default)\n' ...
                       '  2 - AUD\n' ...
                       '  3 - ECN\n' ...
                       '(press Enter to use default [1]): ']);
    if isempty(com_ID); com_ID = 1; end    
    
    %-- Step 2: select ICA order to visualize ---%
    fprintf('Available model orders: %s\n', num2str(modelOrders));
    orderIndex = input(['Select which model order to view'  ...
                        ', press Enter for default [' num2str(modelOrders(end)) ']): ']);
    if isempty(orderIndex); orderIndex = modelOrders(end); end

    %-- Step 3: Display spatial map and change curve ---%
    if analysisMode == 2
        ssm = squeeze(abs(SM_best(find(modelOrders == orderIndex), com_ID, :)));
        noise_ind = find(ssm < 0.5); 
        ssm(noise_ind) = 0;
        ssp = squeeze(angle(SM_best(find(modelOrders == orderIndex), com_ID, :))); 
        ssp(noise_ind) = 0;
        
        ssp_3D = zeros(53,63,46); 
        ssp_3D(inbrain_ind) = ssp;
        ssm_3D = zeros(53,63,46); 
        ssm_3D(inbrain_ind) = ssm;
        figure;SM_display_multiSlice(ssm_3D,0.5,com_ID);

         subRegionOptions = {
            'DMN: 1 - ACC, 2 - PCC, 3 - L-IPL, 4 - R-IPL';
            'AUD: 1 - L-AUD, 2 - R-AUD';
            'ECN: 1 - L-ECN, 2 - M-ECN, 3 - R-ECN'
        };
        fprintf('Available sub-regions for selected network:\n%s\n', subRegionOptions{com_ID});

        defaultSubRegion = 1;
        subRegionIndex = input(['Select one sub-region within the component ' ...
                                '(press Enter for default [' num2str(defaultSubRegion) ']): ']);
        if isempty(subRegionIndex); sr = defaultSubRegion;
        else; sr = subRegionIndex; end
        
        ssm_3D_r = subregion_division(ssm_3D, com_ID, 1);
        region_ssm = squeeze(ssm_3D_r(:,:,:,sr));
        ssp_3D_r = subregion_division(ssp_3D, com_ID, 1);
        region_ssp = squeeze(ssp_3D_r(:,:,:,sr));
        
        [peak_val, idx] = max(region_ssm(:));
        [loc_x, loc_y, loc_z] = ind2sub(size(region_ssm), idx);
        MNI_location = [85 + (loc_x - 25) * 3, 177 + (loc_y - 55) * 3, 63 + (loc_z - 15) * 3];
        
        SMdisplay(region_ssm(inbrain_ind), inbrain_ind, MNI_location);
        SMdisplay(region_ssp(inbrain_ind), inbrain_ind, MNI_location);
        disp_PeakCurve(region_ssm, region_ssp, loc_x, loc_y, loc_z, com_ID, sr);
    else
        mag = squeeze(abs(SM_best(find(modelOrders == orderIndex), com_ID, :)));
        mag_3D = zeros(53,63,46);mag_3D(inbrain_ind) = mag;
        figure;SM_display_multiSlice(mag_3D,1.5,com_ID);

        subRegionOptions = {
        'DMN: 1 - ACC, 2 - PCC, 3 - L-IPL, 4 - R-IPL';
        'AUD: 1 - L-AUD, 2 - R-AUD';
        'ECN: 1 - L-ECN, 2 - M-ECN, 3 - R-ECN'
         };
        fprintf('Available sub-regions for selected network:\n%s\n', subRegionOptions{com_ID});

        defaultSubRegion = 1;
        subRegionIndex = input(['Select one sub-region within the component ' ...
                            '(press Enter for default [' num2str(defaultSubRegion) ']): ']);
        if isempty(subRegionIndex); sr = defaultSubRegion;
        else; sr = subRegionIndex; end
    
        mag_3D_r = subregion_division(mag_3D, com_ID, 1);
        region_mag = squeeze(mag_3D_r(:,:,:,sr));
        
        [peak_val, idx] = max(region_mag(:));
        [loc_x, loc_y, loc_z] = ind2sub(size(region_mag), idx);
        MNI_location = [85 + (loc_x - 25) * 3, 177 + (loc_y - 55) * 3, 63 + (loc_z - 15) * 3];
        
        SMdisplay(region_mag(inbrain_ind), inbrain_ind, MNI_location);
        fprintf('Magnitude-only data: Phase display is not available.\n');
    end

    %% ----------------------- Fitting modeling ------------------------ %%
    w = 5;  % cube size
    d = 2;  % step

    if analysisMode == 2
        SM_concatenate = reshape(squeeze(SM_select(:,:,com_ID,:)), ...
            [size(SM_select,1)*size(SM_select,2), size(SM_select,4)]);

        n_txyz = []; 
        am = [];
        for ii = 1:size(SM_concatenate,1)
            ssm = abs(SM_concatenate(ii,:));
            ssm_3D = zeros(53,63,46); 
            ssm_3D(inbrain_ind) = ssm;
            region_ssm = subregion_division(ssm_3D, com_ID, 0);
            region_ssm = squeeze(region_ssm{sr});

            ssp = angle(SM_concatenate(ii,:));
            ssp_3D = zeros(53,63,46); 
            ssp_3D(inbrain_ind) = ssp;
            region_ssp = subregion_division(ssp_3D, com_ID, 0);
            region_ssp = squeeze(region_ssp{sr});

            [nt_3direction, am(ii,:)] = cal_cube_phase_magnitude(w, d, region_ssp, region_ssm);
            n_txyz(ii,:) = mean(nt_3direction, 1);
        end
        
        [n_txyz_c, a_m_c] = remove_outliers(mean(n_txyz,1).', mean(am,1).');
        
        ShowName = sprintf('%s_ID%03d_%s_%s_mo%s_%s', ...
            personName, infoID, dataTypeStr, icaTypeStr, modelOrdersStr, runStr);
        
        [Params.b0, Params.b1, Params.b2, Params.a_b, Params.n_b, error] = ...
            piecewise_linear_fit(a_m_c, n_txyz_c, 1, 1, ShowName);
    else
    end

    %% -------- Exit / Continue menu --------
    userChoice = input(['\nModeling finished. Next step:\n' ...
                        '  1 - View next region\n' ...
                        '  0 - Exit\n' ...
                        'Enter 1 or 0: ']);

    while ~ismember(userChoice, [0, 1])
        userChoice = input('Invalid selection. Enter 1 (next) or 0 (exit): ');
    end

    if userChoice == 0
        continueFlag = false;
        disp('Exit.');
    end

end


