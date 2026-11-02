function [n_txyz_c, a_m_c] = remove_outliers(n_txyz, a_m)
% remove_outliers - Remove outliers based on interquartile range
%
% Inputs:
%   - n_txyz: mean number of phase polarity transitions along x, y, z directions [1 x N]
%             N = number of phase voxel cubes
%   - a_m: mean activation magnitude [1 x N]
%
% Outputs:
%   - n_txyz_c: n_txyz after removing outliers
%   - a_m_c: a_m after removing outliers

%% identify outliers for n_txyz    
Q1_ntxyz = prctile(n_txyz, 25);  % 25th percentile
Q3_ntxyz = prctile(n_txyz, 75);  % 75th percentile
IQR_ntxyz = Q3_ntxyz - Q1_ntxyz;  % IQR
ol_ntxyz = find((n_txyz > Q3_ntxyz + 1.5 * IQR_ntxyz) | (n_txyz < Q1_ntxyz - 1.5 * IQR_ntxyz));  % outliers

% identify outliers for a_m
Q1_am = prctile(a_m, 25);
Q3_am = prctile(a_m, 75);
IQR_am = Q3_am - Q1_am;
ol_am = find((a_m > (Q3_am + 1.5 * IQR_am)) | (a_m < (Q1_am - 1.5 * IQR_am)));

% Remove outliers
n_txyz_c = n_txyz;
n_txyz_c(union(ol_ntxyz, ol_am)) = [];

a_m_c = a_m;
a_m_c(union(ol_ntxyz, ol_am)) = [];

% Check results
if (max(n_txyz_c) == 0) || (max(a_m_c) == 0)
    indc0 = find(a_m == 0);
    a_m(indc0) = [];
    n_txyz(indc0) = [];
    
    Q1_ntxyz = prctile(n_txyz, 25);
    Q3_ntxyz = prctile(n_txyz, 75);
    IQR_ntxyz = Q3_ntxyz - Q1_ntxyz;
    ol_ntxyz = find((n_txyz > Q3_ntxyz + 1.5 * IQR_ntxyz) | (n_txyz < Q1_ntxyz - 1.5 * IQR_ntxyz));  % outliers

    Q1_am = prctile(a_m, 25);
    Q3_am = prctile(a_m, 75);
    IQR_am = Q3_am - Q1_am;
    ol_am = find((a_m > (Q3_am + 1.5 * IQR_am)) | (a_m < (Q1_am - 1.5 * IQR_am)));

    n_txyz_c = n_txyz;
    n_txyz_c(union(ol_ntxyz, ol_am)) = [];
    a_m_c = a_m;
    a_m_c(union(ol_ntxyz, ol_am)) = [];
end
end
