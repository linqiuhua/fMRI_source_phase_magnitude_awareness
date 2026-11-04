function [performance_mag, performance_pol, performance_fusion] = ...
    SVM_classification(am_multi_order, ntxyz_multi_order)

% SVM_classification - Perform SVM-based classification using voxel cube
% based activation amplitude, phase polarity transition, and their fusion.
%
% Inputs:
%   am_multi_order: Activation amplitude of all voxel cubes
%                   [K x M x N]
%   ntxyz_multi_order: Number of phase polarity transitions
%                      [K x M x N]
%
% Outputs:
%   performance_mag: Classification performance for magnitude features
%   performance_pol: Classification performance for phase features
%   performance_fusion: Classification performance for fused features


%% Sort voxel cubes by activation amplitude (descending) per subject

for sub_ind = 1:size(am_multi_order,1)

    [~, am_order] = sort( ...
        squeeze(mean(am_multi_order(sub_ind,:,:), 2)), 'descend');

    ordered_amplitude(sub_ind,:,:) = ...
        squeeze(am_multi_order(sub_ind,:,am_order));

    ordered_polarity(sub_ind,:,:) = ...
        squeeze(ntxyz_multi_order(sub_ind,:,am_order));

end


%% Classification using activation amplitude

X_hc = reshape(ordered_amplitude(1:40,:,:), ...
    [40*size(ordered_amplitude,2), size(ordered_amplitude,3)]);

X_sz = reshape(ordered_amplitude(41:82,:,:), ...
    [42*size(ordered_amplitude,2), size(ordered_amplitude,3)]);

X_all = [X_hc; X_sz];

y_all = [ones(40*size(ordered_amplitude,2),1); ...
         zeros(42*size(ordered_amplitude,2),1)];


% Five-fold outer cross-validation
cv_fold = cvpartition(y_all, 'KFold', 5);

box_vals = logspace(-3, 3, 7);
scale_vals = logspace(-3, 3, 7);

performance_mag = struct( ...
    'accuracy', [], 'sensitivity', [], 'specificity', []);

for fold = 1:cv_fold.NumTestSets

    % Outer training and testing sets
    train_idx = training(cv_fold, fold);
    test_idx = test(cv_fold, fold);

    X_train = X_all(train_idx,:);
    y_train = y_all(train_idx);

    X_test = X_all(test_idx,:);
    y_test = y_all(test_idx);


    %% Hyperparameter optimization
    cv_inner = cvpartition(y_train, 'KFold', 5);

    best_cv_accuracy = 0;
    for i = 1:length(box_vals)
        for j = 1:length(scale_vals)
            C_val = box_vals(i);
            scale_val = scale_vals(j);
            acc = zeros(cv_inner.NumTestSets,1);
            for inner = 1:cv_inner.NumTestSets % 20% validation
                inner_train_idx = training(cv_inner, inner);
                inner_val_idx = test(cv_inner, inner);
                model = fitcsvm( ... % self define
                    X_train(inner_train_idx,:), ... 
                    y_train(inner_train_idx), ...
                    'KernelFunction', 'rbf', ... 
                    'BoxConstraint', C_val, ...
                    'KernelScale', scale_val, ...
                    'Standardize', false);

                y_pred = predict(model, ...
                    X_train(inner_val_idx,:));

                acc(inner) = mean( ...
                    y_pred == y_train(inner_val_idx));

            end

            mean_acc = mean(acc);

            if mean_acc > best_cv_accuracy

                best_cv_accuracy = mean_acc;

                best_params.C = C_val;
                best_params.kernelScale = scale_val;

            end

        end
    end
    model = fitcsvm(X_train, y_train, ...
        'KernelFunction', 'rbf', ...
        'BoxConstraint', best_params.C, ...
        'KernelScale', best_params.kernelScale, ...
        'Standardize', false);

    y_pred = predict(model, X_test); % 20% test

    TP = sum((y_test == 1) & (y_pred == 1));
    TN = sum((y_test == 0) & (y_pred == 0));
    FP = sum((y_test == 0) & (y_pred == 1));
    FN = sum((y_test == 1) & (y_pred == 0));

    performance_mag.accuracy(fold) = ...
        (TP + TN) / (TP + TN + FP + FN) * 100;

    performance_mag.sensitivity(fold) = ...
        TP / (TP + FN) * 100;

    performance_mag.specificity(fold) = ...
        TN / (TN + FP) * 100;

end


%% Classification using phase polarity transitions

X_hc = reshape(ordered_polarity(1:40,:,:), ...
    [40*size(ordered_amplitude,2), size(ordered_polarity,3)]);

X_sz = reshape(ordered_polarity(41:82,:,:), ...
    [42*size(ordered_amplitude,2), size(ordered_polarity,3)]);

X_all = [X_hc; X_sz];

y_all = [ones(40*size(ordered_amplitude,2),1); ...
         zeros(42*size(ordered_amplitude,2),1)];


% Five-fold outer cross-validation
cv_fold = cvpartition(y_all, 'KFold', 5);

performance_pol = struct( ...
    'accuracy', [], 'sensitivity', [], 'specificity', []);

for fold = 1:cv_fold.NumTestSets

    train_idx = training(cv_fold, fold);
    test_idx = test(cv_fold, fold);

    X_train = X_all(train_idx,:);
    y_train = y_all(train_idx);

    X_test = X_all(test_idx,:);
    y_test = y_all(test_idx);

    cv_inner = cvpartition(y_train, 'KFold', 5);

    best_cv_accuracy = 0;

    for i = 1:length(box_vals)
        for j = 1:length(scale_vals)

            C_val = box_vals(i);
            scale_val = scale_vals(j);

            acc = zeros(cv_inner.NumTestSets,1);

            for inner = 1:cv_inner.NumTestSets

                inner_train_idx = training(cv_inner, inner);
                inner_val_idx = test(cv_inner, inner);

                model = fitcsvm( ...
                    X_train(inner_train_idx,:), ...
                    y_train(inner_train_idx), ...
                    'KernelFunction', 'rbf', ...
                    'BoxConstraint', C_val, ...
                    'KernelScale', scale_val, ...
                    'Standardize', false);

                y_pred = predict(model, ...
                    X_train(inner_val_idx,:));

                acc(inner) = mean( ...
                    y_pred == y_train(inner_val_idx));

            end

            mean_acc = mean(acc);

            if mean_acc > best_cv_accuracy

                best_cv_accuracy = mean_acc;

                best_params.C = C_val;
                best_params.kernelScale = scale_val;

            end

        end
    end

    model = fitcsvm(X_train, y_train, ...
        'KernelFunction', 'rbf', ...
        'BoxConstraint', best_params.C, ...
        'KernelScale', best_params.kernelScale, ...
        'Standardize', false);

    y_pred = predict(model, X_test);

    TP = sum((y_test == 1) & (y_pred == 1));
    TN = sum((y_test == 0) & (y_pred == 0));
    FP = sum((y_test == 0) & (y_pred == 1));
    FN = sum((y_test == 1) & (y_pred == 0));

    performance_pol.accuracy(fold) = ...
        (TP + TN) / (TP + TN + FP + FN) * 100;

    performance_pol.sensitivity(fold) = ...
        TP / (TP + FN) * 100;

    performance_pol.specificity(fold) = ...
        TN / (TN + FP) * 100;

end


%% Classification using fusion features

fusion_feature_matrix = cat(3, ordered_amplitude, ordered_polarity);

X_hc = reshape(fusion_feature_matrix(1:40,:,:), ...
    [40*size(ordered_amplitude,2), size(fusion_feature_matrix,3)]);

X_sz = reshape(fusion_feature_matrix(41:82,:,:), ...
    [42*size(ordered_amplitude,2), size(fusion_feature_matrix,3)]);

X_all = [X_hc; X_sz];

y_all = [ones(40*size(ordered_amplitude,2),1); ...
         zeros(42*size(ordered_amplitude,2),1)];


% Five-fold outer cross-validation
cv_fold = cvpartition(y_all, 'KFold', 5);

performance_fusion = struct( ...
    'accuracy', [], 'sensitivity', [], 'specificity', []);

for fold = 1:cv_fold.NumTestSets

    train_idx = training(cv_fold, fold);
    test_idx = test(cv_fold, fold);

    X_train = X_all(train_idx,:);
    y_train = y_all(train_idx);

    X_test = X_all(test_idx,:);
    y_test = y_all(test_idx);

    cv_inner = cvpartition(y_train, 'KFold', 5);

    best_cv_accuracy = 0;

    for i = 1:length(box_vals)
        for j = 1:length(scale_vals)

            C_val = box_vals(i);
            scale_val = scale_vals(j);

            acc = zeros(cv_inner.NumTestSets,1);

            for inner = 1:cv_inner.NumTestSets

                inner_train_idx = training(cv_inner, inner);
                inner_val_idx = test(cv_inner, inner);

                model = fitcsvm( ...
                    X_train(inner_train_idx,:), ...
                    y_train(inner_train_idx), ...
                    'KernelFunction', 'rbf', ...
                    'BoxConstraint', C_val, ...
                    'KernelScale', scale_val, ...
                    'Standardize', false);

                y_pred = predict(model, ...
                    X_train(inner_val_idx,:));

                acc(inner) = mean( ...
                    y_pred == y_train(inner_val_idx));

            end

            mean_acc = mean(acc);

            if mean_acc > best_cv_accuracy

                best_cv_accuracy = mean_acc;

                best_params.C = C_val;
                best_params.kernelScale = scale_val;

            end

        end
    end

    model = fitcsvm(X_train, y_train, ...
        'KernelFunction', 'rbf', ...
        'BoxConstraint', best_params.C, ...
        'KernelScale', best_params.kernelScale, ...
        'Standardize', false);

    y_pred = predict(model, X_test);

    TP = sum((y_test == 1) & (y_pred == 1));
    TN = sum((y_test == 0) & (y_pred == 0));
    FP = sum((y_test == 0) & (y_pred == 1));
    FN = sum((y_test == 1) & (y_pred == 0));

    performance_fusion.accuracy(fold) = ...
        (TP + TN) / (TP + TN + FP + FN) * 100;

    performance_fusion.sensitivity(fold) = ...
        TP / (TP + FN) * 100;

    performance_fusion.specificity(fold) = ...
        TN / (TN + FP) * 100;

end

end
