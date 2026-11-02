%**************************************************************************
% Description:
%   Group-level analysis code for "Uncovering abnormal patterns of schizophrenia from fMRI regional source phase"
% submitted to Nature Computational Science%
%
% Functionality:
%   - Load complex-valued spatial maps from healthy controls (HC) and schizophrenia (SZ) groups.
%   - Compute activation magnitude and phase polarity transitions within brain sub-region.
%   - Perform piecewise linear fitting to characterize the relationship
%     between phase polarity transitions and activation magnitude.
%   - Compare the five model parameters and fitting errors between HC and SZ
%     groups using two-sample t-tests.
%   - Evaluate group differences with BH-FDR correction.
%   - Identify and visualize biomarker region with significant group difference and their spatial overlap.
%   - Evaluate the discriminative performance of individual features and their fusion using SVM.
%
% Input:
%   - sm_COBRE82subjects_DMN.mat for example:
%       Complex-valued spatial maps of the DMN component from 82 subjects
%       in the COBRE dataset.
%       Sm: complex-valued spatial maps [M x K x V]
%              M = number of model orders
%              K = number of subjects
%              V = number of in-brain voxels
%
%   - inbrain_ind.mat: indices of in-brain voxels [1 x V]
%              V = number of in-brain voxels
%
%   - Analysis settings:
%       w       : Cube size.
%       d       : Sliding step.
%       com_ID  : Index of the selected functional network. 1 - DMN (Default Mode Network)
%       sr      : Index of the selected brain sub-region. 1 - ACC 
%
% Output:
%   - Params:
%       Structure containing the fitted piecewise linear model parameters:
%           .b0 : Intercept
%           .b1 : Slope of the first segment
%           .b2 : Slope of the second segment
%           .a_b: Breakpoint of activation magnitude
%           .n_b: Breakpoint of phase polarity transitions
%
%   - error:
%       Fitting residuals (MSE) of the piecewise linear model.
%
%   - FDR_result:
%       Structure containing the original and BH-FDR-corrected
%       significant cube indices and statistical results for both features.
%
%   - Spatial visualization:
%       Significant activation magnitude regions, phase polarity regions,
%       and their shared and distinct spatial maps.
%
%   - SVM classification results:
%       Accuracy, sensitivity, and specificity of activation magnitude,
%       phase polarity transitions, and their fusion.
%
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

%% ---------- Load complex-valued spatial maps for the HC and SZ groups ---------- %%
fprintf('Loading complex-valued spatial maps for the DMN component (or optionally using our precomputed results)...\n');
load data\sm_COBRE82subjects_DMN.mat
load data\inbrain_ind.mat;
%% --------------------------Individual fitting modeling ------------------------- %%
w = 5;  % cube size
d = 2;  % step

com_ID = 1; % DMN component
sr = 1; % ACC sub-region

for subj=1:size(Sm,2)
    fprintf('Analyzing subject %d/%d...\n', subj, size(Sm, 2));
    n_txyz = []; am = [];
    for ii = 1:size(Sm,1)
        ssm = squeeze(abs(Sm(ii,subj,:)));
        noise_ind = find(ssm < 0.5); 
        ssm(noise_ind) = 0;
        ssm_3D = zeros(53,63,46); 
        ssm_3D(inbrain_ind) = ssm;
        region_ssm = subregion_division(ssm_3D, com_ID, 0);
        region_ssm = squeeze(region_ssm{sr});
    
        ssp = squeeze(angle(Sm(ii,subj,:)));
        ssp(noise_ind) = 0;
        ssp_3D = zeros(53,63,46); 
        ssp_3D(inbrain_ind) = ssp;
        region_ssp = subregion_division(ssp_3D, com_ID, 0);
        region_ssp = squeeze(region_ssp{sr});
    
        [nt_3direction, am(ii,:)] = cal_cube_phase_magnitude(w, d, region_ssp, region_ssm);
        n_txyz(ii,:) = mean(nt_3direction, 1);
    end
    n_txyz_all(subj,:,:) = n_txyz;am_all(subj,:,:) = am;
    [n_txyz_c, a_m_c] = remove_outliers(mean(n_txyz,1).', mean(am,1).');
    
    [Params.b0(subj), Params.b1(subj), Params.b2(subj), Params.a_b(subj), Params.n_b(subj), error(subj)] = ...
        piecewise_linear_fit(a_m_c, n_txyz_c, 0, 1, '');
end

%% ------------------------ Group statistical analyses ------------------------- %%
% Set group labels for healthy controls and schizophrenia patients
if size(Sm,2) == 82 % COBRE dataset
    HC_ind = 1:40;       % healthy controls
    SZ_ind = 41:82;      % schizophrenia patients
elseif size(Sm,2) == 100 % UCLA dataset
    HC_ind = 1:50;       % healthy controls
    SZ_ind = 51:100;     % schizophrenia patients

else
    error('Please set the subject labels for the HC and SZ groups.');
end

%--step 1: group difference in piecewise linear model parameters ----%
[~,pvalue_b0,~,stat]=ttest2(Params.b0(SZ_ind),Params.b0(HC_ind));
figure;subplot(2,3,1);Boxplot(Params.b0(SZ_ind).', Params.b0(HC_ind).',pvalue_b0);title('b0')
[~,pvalue_b1,~,stat]=ttest2(Params.b1(SZ_ind),Params.b1(HC_ind));
subplot(2,3,2);Boxplot(Params.b1(SZ_ind).', Params.b1(HC_ind).',pvalue_b1);title('b1')
[~,pvalue_b2,~,stat]=ttest2(Params.b2(SZ_ind),Params.b2(HC_ind));
subplot(2,3,3);Boxplot(Params.b2(SZ_ind).', Params.b2(HC_ind).',pvalue_b2);title('b2')
[~,pvalue_ab,~,stat]=ttest2(Params.a_b(SZ_ind),Params.a_b(HC_ind));
subplot(2,3,4);Boxplot(Params.a_b(SZ_ind).', Params.a_b(HC_ind).',pvalue_ab);title('ab')
[~,pvalue_nb,~,stat]=ttest2(Params.n_b(SZ_ind),Params.n_b(HC_ind));
subplot(2,3,5);Boxplot(Params.n_b(SZ_ind).', Params.n_b(HC_ind).',pvalue_nb);title('nb')
% fitting errors
subplot(2,3,6);Boxplot(error(SZ_ind).', error(HC_ind).',Inf);ylim([0 0.06]);title('MSE')
% Overall figure title
sgt = sgtitle('Group differences in five model parameters', ...
    'FontSize', 16, 'FontWeight', 'normal');
% Adjust the vertical position of the top-row subplots
ax = findall(gcf, 'Type', 'axes');
for k = 1:length(ax)
    pos = get(ax(k), 'Position');

    if pos(2) > 0.5
        pos(2) = pos(2) - 0.04;
        set(ax(k), 'Position', pos);
    end
end

%--step 2: group difference in phase polarity transitionse and activation magnitude ----%
% rank magnitude and phase cubes 
n_txyz_mean = squeeze(mean(n_txyz_all,2)); % average of multiple model orders
am_mean = squeeze(mean(am_all,2));
for subj=1:size(n_txyz_mean,1)
    [order_am(subj,:),order_ind]=sort(am_mean(subj,:),'descend');
    order_ntxyz(subj,:)=n_txyz_mean(subj,order_ind);
    common_ind(subj,:)=numel(find(order_am(subj,:)==0));   
end
% evaluate group differences at each voxel index
for cube_ind=1:size(order_am,2)
    [~,pvalue_am(cube_ind),~,stat_am]=ttest2(order_am(SZ_ind,cube_ind),order_am(HC_ind,cube_ind));
    tvalue_am(cube_ind)=stat_am.tstat;
    [~,pvalue_ntxyz(cube_ind),~,stat_ntxyz]=ttest2(order_ntxyz(SZ_ind,cube_ind),order_ntxyz(HC_ind,cube_ind));
    tvalue_ntxyz(cube_ind)=stat_ntxyz.tstat;
end
FDR_result = FDR_group_differences( pvalue_am,pvalue_ntxyz, order_am, common_ind, ...
      0.05,1,length(tvalue_am), tvalue_am, tvalue_ntxyz);
figure;subplot(2,1,1);
plot_tvalue=zeros(size(tvalue_am));
plot_tvalue(:,find(pvalue_am>=0.05))=tvalue_am(:,find(pvalue_am>=0.05));
bar(plot_tvalue,1,'Facecolor',[222,242,237]/255,'Edgecolor',[222,242,237]/255);
plot_tvalue=zeros(size(tvalue_am));plot_tvalue(:,find(pvalue_am<0.05))=tvalue_am(:,find(pvalue_am<0.05));
hold on;bar(plot_tvalue,1,'Facecolor',[96,189,170]/255,'Edgecolor',[96,189,170]/255);
hold on; line([(size(tvalue_am,2)-max(common_ind)),(size(tvalue_am,2)-max(common_ind))],[-6,6],'LineStyle',':','LineWidth',1.5,'Color',[144,144,144]/255);
set(gca,'FontName','Arial','FontSize',20,'LineWidth',0.8, 'XColor',[0 0 0], 'YColor',[0 0 0], 'Box','on');ylim([-6 6])
plot_tvalue=zeros(size(tvalue_am));plot_tvalue(:,FDR_result.sig_am_fdr)=tvalue_am(:,FDR_result.sig_am_fdr);
hold on;bar(plot_tvalue,1,'Facecolor',[255,175, 202] / 255,'Edgecolor',[255,175, 202] / 255);
line([0,ceil(size(n_txyz,2)/100)*100],[0,0],'LineWidth',0.7,'Color',[100,100,100]/255);xlim([0,ceil(size(n_txyz,2)/100)*100]);
subplot(2,1,2);plot_tvalue=zeros(size(tvalue_ntxyz));plot_tvalue(:,find(pvalue_ntxyz>=0.05))=tvalue_ntxyz(:,find(pvalue_ntxyz>=0.05));
bar(plot_tvalue,1,'Facecolor',[233,230,238]/255,'Edgecolor',[233,230,238]/255);
plot_tvalue=zeros(size(tvalue_ntxyz));
plot_tvalue(:,find(pvalue_ntxyz<0.05))=tvalue_ntxyz(:,find(pvalue_ntxyz<0.05));
hold on;bar(plot_tvalue,1,'Facecolor',[185,170,202]/255,'Edgecolor',[185,170,202]/255);
hold on; line([(size(tvalue_am,2)-max(common_ind)),(size(tvalue_am,2)-max(common_ind))],[-6,6],'LineStyle',':','LineWidth',1.5,'Color',[144,144,144]/255);
set(gca,'FontName','Arial','FontSize',20,'LineWidth',0.8, 'XColor',[0 0 0], 'YColor',[0 0 0], 'Box','on');ylim([-6 6])
plot_tvalue=zeros(size(tvalue_ntxyz));plot_tvalue(:,FDR_result.sig_ntxyz_fdr)=tvalue_ntxyz(:,FDR_result.sig_ntxyz_fdr);
hold on;bar(plot_tvalue,1,'Facecolor',[255, 153, 51] / 255,'Edgecolor',[255, 153, 51] / 255);
line([0,ceil(size(n_txyz,2)/100)*100],[0,0],'LineWidth',0.7,'Color',[100,100,100]/255);xlim([0,ceil(size(n_txyz,2)/100)*100]);


%--step 3: identify the phase and magnitude cubes exhibiting significant group differences ----%
ssm=abs(Sm(10,1,:)); % exampling HC subject
noise_ind=find(ssm<0.5);ssm(noise_ind)=0;
ssm_3D=zeros(53,63,46);ssm_3D(inbrain_ind)=ssm;
ssm_3D_r=subregion_division(ssm_3D,com_ID,0);
region_ssm=squeeze(ssm_3D_r{sr});
ssp=angle(Sm(10,1,:));ssp(noise_ind)=0;
ssp_3D=zeros(53,63,46);ssp_3D(inbrain_ind)=ssp;
ssp_3D_r=subregion_division(ssp_3D,com_ID,0);
region_ssp=squeeze(ssp_3D_r{sr});

% identify significant SZ-HC group difference (without FDR correction)
ntxyz_sig_ind=FDR_result.sig_ntxyz_uncorrected;
ntxyz_sig_tv=tvalue_ntxyz(ntxyz_sig_ind);
am_sig_ind=FDR_result.sig_am_uncorrected;
am_sig_tv=tvalue_am(am_sig_ind);
[sig_ssm,sig_ssp]=...
   sig_cube_phase_magnitude(w,d,region_ssp,region_ssm,am_sig_ind,ntxyz_sig_ind,am_sig_tv,ntxyz_sig_tv);
% mask based on activation range
x_idx = 16:36; y_idx = 43:63; z_idx = 13:30; % ACC sub-region index
region_ssm_3D=zeros(53,63,46);region_ssm_3D(x_idx,y_idx,z_idx)=region_ssm;
ind0=(find(region_ssm_3D==0));
sig_ssm_3D=zeros(53,63,46);sig_ssm_3D(x_idx,y_idx,z_idx)=sig_ssm;sig_ssm_3D(ind0)=0;
sig_ssp_3D=zeros(53,63,46);sig_ssp_3D(x_idx,y_idx,z_idx)=sig_ssp;sig_ssp_3D(ind0)=0;
% calculate shared and distinct region
ssm_sig_ind=find(sig_ssm_3D);ssp_sig_ind=find(sig_ssp_3D);
shared_ind=intersect(ssm_sig_ind,ssp_sig_ind);
unique_am_op = setdiff(ssm_sig_ind, ssp_sig_ind);unique_ntxyz_op = setdiff(ssp_sig_ind, ssm_sig_ind);
over_map=zeros(53,63,46);over_map(shared_ind)=10^(-4);
over_map(unique_am_op)=1;over_map(unique_ntxyz_op)=-1;
over_map(1,10,1)=1;over_map(1,12,1)=-1;

% identify significant SZ-HC group difference (with FDR correction)
ntxyz_sig_ind=FDR_result.sig_ntxyz_fdr;ntxyz_sig_tv=tvalue_ntxyz(ntxyz_sig_ind);
am_sig_ind=FDR_result.sig_am_fdr;am_sig_tv=tvalue_am(am_sig_ind);
[sig_ssm_fdr,sig_ssp_fdr]=...
   sig_cube_phase_magnitude(w,d,region_ssp,region_ssm,am_sig_ind,ntxyz_sig_ind,am_sig_tv,ntxyz_sig_tv);
% masking based on activation range
region_ssm_3D=zeros(53,63,46);region_ssm_3D(x_idx,y_idx,z_idx)=region_ssm;
ind0=(find(region_ssm_3D==0));
sig_ssm_3D=zeros(53,63,46);sig_ssm_3D(x_idx,y_idx,z_idx)=sig_ssm_fdr;sig_ssm_3D(ind0)=0;
sig_ssp_3D=zeros(53,63,46);sig_ssp_3D(x_idx,y_idx,z_idx)=sig_ssp_fdr;sig_ssp_3D(ind0)=0;
% calculate shared and distinct region
ssm_sig_ind_fdr=find(sig_ssm_3D);ssp_sig_ind_fdr=find(sig_ssp_3D);
shared_ind_fdr=intersect(ssm_sig_ind_fdr,ssp_sig_ind_fdr);
unique_am_op_fdr = setdiff(ssm_sig_ind_fdr, ssp_sig_ind_fdr);
unique_ntxyz_op_fdr = setdiff(ssp_sig_ind_fdr, ssm_sig_ind_fdr);
% display
over_lap_all=over_map;
over_lap_all(shared_ind_fdr)=0.8;
over_lap_all(unique_am_op_fdr)=0.4;
over_lap_all(unique_ntxyz_op_fdr)=-0.1;
location = [83,177,65] ;
SMdisplay((sig_ssm_3D(inbrain_ind)), inbrain_ind,location); % biomarker regions formed by the source magnitude
SMdisplay((sig_ssp_3D(inbrain_ind)), inbrain_ind,location); % biomarker regions formed by the source phase polarity
SMdisplay_overlap_fdr((over_lap_all(inbrain_ind)), inbrain_ind, location);% overlap

%--step 4: SVM classification ----%
[measure_am,measure_ntxyz,measure_fusion] = SVM_classification (am_all,n_txyz_all);

% print results
measure_names = {'Activation mangitude', 'Phase Polarity Transitions', 'Fusion'};
measure_structs = {measure_am, measure_ntxyz, measure_fusion};
for i = 1:3
    fprintf('%s (five-fold):\n', measure_names{i});
    fprintf('  Accuracy   : %.2f%% ± %.2f%%\n', mean(measure_structs{i}.accuracy), std(measure_structs{i}.accuracy));
    fprintf('  Sensitivity: %.2f%% ± %.2f%%\n', mean(measure_structs{i}.sensitivity), std(measure_structs{i}.sensitivity));
    fprintf('  Specificity: %.2f%% ± %.2f%%\n\n', mean(measure_structs{i}.specificity), std(measure_structs{i}.specificity));
end



%% ------------------------ Group-level Tmap ------------------------- %%
HC_data = squeeze(abs(Sm(5,HC_ind, :)));  
SZ_data = squeeze(abs(Sm(5,SZ_ind, :))); 
HC_data(HC_data<0.5)=0;
SZ_data(SZ_data<0.5)=0;
n_voxels = size(Sm, 3);

% perform voxel-wise one sample t-test
p_HC = zeros(1, n_voxels);p_SZ = zeros(1, n_voxels);
t_HC = zeros(1, n_voxels);t_SZ = zeros(1, n_voxels);

for vox = 1:n_voxels
    [h, p_HC(vox), ci, stats] = ttest(HC_data(:, vox));
    t_HC(vox) = stats.tstat;
    [h, p_SZ(vox), ci, stats] = ttest(SZ_data(:, vox));
    t_SZ(vox) = stats.tstat;
end
t_HC(p_HC>=0.001)=0;t_SZ(p_SZ>=0.001)=0;

% ----- FDR correction ----- %
q_HC = mafdr(p_HC, 'BHFDR', true);
q_SZ = mafdr(p_SZ, 'BHFDR', true);

alpha_fdr = 0.001;
significant_HC = q_HC < alpha_fdr;
significant_SZ = q_SZ < alpha_fdr;

HC_sig_idx = find(significant_HC);
SZ_sig_idx = find(significant_SZ);
fdr_hc=t_HC.*significant_HC;
fdr_sz=t_SZ.*significant_SZ;

sm_hc=zeros(53,63,46); sm_hc(inbrain_ind)=fdr_hc;
sm_hc(isnan(sm_hc))=0;

sm_sz=zeros(53,63,46); sm_sz(inbrain_ind)=fdr_sz;
sm_sz(isnan(sm_sz))=0;

figure;SM_display_multiSlice(sm_hc,0,1);
figure;SM_display_multiSlice(sm_sz,0,1);
