function [performance_mag, performance_pol, performance_fusion] = SVM_classification(am_multi_order, ntxyz_multi_order)
% SVM_classification - Perform SVM-based classification using voxel cube
% based activation amplitude, phase polarity transition, and their fusion.
%
% Inputs:
%   - am_multi_order: activation amplitude of all voxel cubes for multiple model orders [K x M x N]
%        K = number of subjects; M = number of model orders; N = number of magnitude voxel cubes
%   - ntxyz_multi_order: number of phase polarity transitions of all voxel cubes
%        for multiple model orders [K x M x N]
%
% Outputs:
%   - performance_mag: struct with fields: accuracy, sensitivity, specificity
%   - performance_pol: same as above, for polarity features
%   - performance_fusion: same as above, for fused features

%% Sort voxel cubes by activation amplitude (descending) per subject
for sub_ind = 1:size(am_multi_order,1)
    [~, am_order] = sort(squeeze(mean(am_multi_order(sub_ind,:,:), 2)), 'descend');
    ordered_amplitude(sub_ind,:,:) = squeeze(am_multi_order(sub_ind,:,am_order));
    ordered_polarity(sub_ind,:,:)  = squeeze(ntxyz_multi_order(sub_ind,:,am_order));
end

%% classification using activation amplitude
X_hc = reshape(ordered_amplitude(1:40,:,:), [40*size(ordered_amplitude,2), size(ordered_amplitude,3)]);
X_sz = reshape(ordered_amplitude(41:82,:,:), [42*size(ordered_amplitude,2), size(ordered_amplitude,3)]);
X_all = [X_hc; X_sz];
y_all = [ones(40*size(ordered_amplitude,2),1); zeros(42*size(ordered_amplitude,2),1)];

% five-fold cross validation
cv_fold = cvpartition(y_all, 'KFold', 5);
box_vals1 = logspace(-3, 3, 7);
scale_vals1 = logspace(-3, 3, 7);
box_vals = logspace(-3, 3, 7);
scale_vals = logspace(-3, 3, 7);
best_cv_accuracy = 0;

for i = 1:length(box_vals1)
    for j = 1:length(scale_vals1)
        C_val = box_vals1(i);
        scale_val = scale_vals1(j);
        acc = zeros(cv_fold.NumTestSets,1);
        for fold = 1:cv_fold.NumTestSets
            train_idx = training(cv_fold, fold);
            test_idx = test(cv_fold, fold);
            model = fitcsvm(X_all(train_idx,:), y_all(train_idx), 'KernelFunction', 'rbf', ...
                'BoxConstraint', C_val, 'KernelScale', scale_val, 'Standardize', false);
            y_pred = predict(model, X_all(test_idx,:));
            acc(fold) = mean(y_pred == y_all(test_idx));
        end
        mean_acc = mean(acc);
        if mean_acc > best_cv_accuracy
            best_cv_accuracy = mean_acc;
            best_params.C = C_val;
            best_params.kernelScale = scale_val;
        end
    end
end

% test using optimal parameters
performance_mag = struct('accuracy', [], 'sensitivity', [], 'specificity', []);
for fold = 1:cv_fold.NumTestSets
    train_idx = training(cv_fold, fold);
    test_idx  = test(cv_fold, fold);
    model = fitcsvm(X_all(train_idx,:), y_all(train_idx), 'KernelFunction', 'rbf', ...
        'BoxConstraint', best_params.C, 'KernelScale', best_params.kernelScale, 'Standardize', false);
    y_test = y_all(test_idx);
    y_pred = predict(model, X_all(test_idx,:));
    TP = sum((y_test == 1) & (y_pred == 1));
    TN = sum((y_test == 0) & (y_pred == 0));
    FP = sum((y_test == 0) & (y_pred == 1));
    FN = sum((y_test == 1) & (y_pred == 0));
    performance_mag.accuracy(fold) = (TP + TN) / (TP + TN + FP + FN) * 100;
    performance_mag.sensitivity(fold) = TP / (TP + FN) * 100;
    performance_mag.specificity(fold) = TN / (TN + FP) * 100;
end

%% classification using phase polarity transitions
X_hc = reshape(ordered_polarity(1:40,:,:), [40*size(ordered_amplitude,2), size(ordered_polarity,3)]);
X_sz = reshape(ordered_polarity(41:82,:,:), [42*size(ordered_amplitude,2), size(ordered_polarity,3)]);
X_all = [X_hc; X_sz];
y_all = [ones(40*size(ordered_amplitude,2),1); zeros(42*size(ordered_amplitude,2),1)];

cv_fold = cvpartition(y_all, 'KFold', 5);
best_cv_accuracy = 0;

for i = 1:length(box_vals)
    for j = 1:length(scale_vals)
        C_val = box_vals(i);
        scale_val = scale_vals(j);
        acc = zeros(cv_fold.NumTestSets,1);
        for fold = 1:cv_fold.NumTestSets
            train_idx = training(cv_fold, fold);
            test_idx = test(cv_fold, fold);
            model = fitcsvm(X_all(train_idx,:), y_all(train_idx), 'KernelFunction', 'rbf', ...
                'BoxConstraint', C_val, 'KernelScale', scale_val, 'Standardize', false);
            y_pred = predict(model, X_all(test_idx,:));
            acc(fold) = mean(y_pred == y_all(test_idx));
        end
        mean_acc = mean(acc);
        if mean_acc > best_cv_accuracy
            best_cv_accuracy = mean_acc;
            best_params.C = C_val;
            best_params.kernelScale = scale_val;
        end
    end
end

performance_pol = struct('accuracy', [], 'sensitivity', [], 'specificity', []);
for fold = 1:cv_fold.NumTestSets
    train_idx = training(cv_fold, fold);
    test_idx  = test(cv_fold, fold);
    model = fitcsvm(X_all(train_idx,:), y_all(train_idx), 'KernelFunction', 'rbf', ...
        'BoxConstraint', best_params.C, 'KernelScale', best_params.kernelScale, 'Standardize', false);
    y_test = y_all(test_idx);
    y_pred = predict(model, X_all(test_idx,:));
    TP = sum((y_test == 1) & (y_pred == 1));
    TN = sum((y_test == 0) & (y_pred == 0));
    FP = sum((y_test == 0) & (y_pred == 1));
    FN = sum((y_test == 1) & (y_pred == 0));
    performance_pol.accuracy(fold) = (TP + TN) / (TP + TN + FP + FN) * 100;
    performance_pol.sensitivity(fold) = TP / (TP + FN) * 100;
    performance_pol.specificity(fold) = TN / (TN + FP) * 100;
end

%% Classification using fusion features
fusion_feature_matrix = cat(3, ordered_amplitude, ordered_polarity);

X_hc = reshape(fusion_feature_matrix(1:40,:,:), [40*size(ordered_amplitude,2), size(fusion_feature_matrix,3)]);
X_sz = reshape(fusion_feature_matrix(41:82,:,:), [42*size(ordered_amplitude,2), size(fusion_feature_matrix,3)]);
X_all = [X_hc; X_sz];
y_all = [ones(40*size(ordered_amplitude,2),1); zeros(42*size(ordered_amplitude,2),1)];

cv_fold = cvpartition(y_all, 'KFold', 5);
best_cv_accuracy = 0;

for i = 1:length(box_vals)
    for j = 1:length(scale_vals)
        C_val = box_vals(i);
        scale_val = scale_vals(j);
        acc = zeros(cv_fold.NumTestSets,1);
        for fold = 1:cv_fold.NumTestSets
            train_idx = training(cv_fold, fold);
            test_idx = test(cv_fold, fold);
            model = fitcsvm(X_all(train_idx,:), y_all(train_idx), 'KernelFunction', 'rbf', ...
                'BoxConstraint', C_val, 'KernelScale', scale_val, 'Standardize', false);
            y_pred = predict(model, X_all(test_idx,:));
            acc(fold) = mean(y_pred == y_all(test_idx));
        end
        mean_acc = mean(acc);
        if mean_acc > best_cv_accuracy
            best_cv_accuracy = mean_acc;
            best_params.C = C_val;
            best_params.kernelScale = scale_val;
        end
    end
end

performance_fusion = struct('accuracy', [], 'sensitivity', [], 'specificity', []);
for fold = 1:cv_fold.NumTestSets
    train_idx = training(cv_fold, fold);
    test_idx  = test(cv_fold, fold);
    model = fitcsvm(X_all(train_idx,:), y_all(train_idx), 'KernelFunction', 'rbf', ...
        'BoxConstraint', best_params.C, 'KernelScale', best_params.kernelScale, 'Standardize', false);
    y_test = y_all(test_idx);
    y_pred = predict(model, X_all(test_idx,:));
    TP = sum((y_test == 1) & (y_pred == 1));
    TN = sum((y_test == 0) & (y_pred == 0));
    FP = sum((y_test == 0) & (y_pred == 1));
    FN = sum((y_test == 1) & (y_pred == 0));
    performance_fusion.accuracy(fold) = (TP + TN) / (TP + TN + FP + FN) * 100;
    performance_fusion.sensitivity(fold) = TP / (TP + FN) * 100;
    performance_fusion.specificity(fold) = TN / (TN + FP) * 100;
end

end
