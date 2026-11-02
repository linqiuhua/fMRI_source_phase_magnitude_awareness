function [SM_interest, SM_best] = selection_component_real (SM, TC, reference)
%selection_component - Select the component of interest based on spatial correlation with reference.
%   Inputs:
%    - SM: cell array of spatial components [Num_run x Num_order],
%          each cell contains [MO x V] real-valued spatial maps.
%    - TC: cell array of time courses [Num_run x Num_order],
%          each cell contains [T x MO] complex-valued time courses.
%    - reference: vector of spatial reference template [1 x V].
%
%   Output:
%    - SM_interest: selected component [Num_run x Num_order x V].
%    - SM_best: best component [Num_order x V].
%% Component selection
for r = 1:size(SM, 1)  % Loop over ICA runs
    for m = 1:size(SM, 2)  % Loop over model orders
        Se = SM{r, m};       % Spatial components [MO x V]
        At = TC{r, m};       % Time courses [T x MO]

        num_comp = size(Se, 1);
        cc = zeros(1, num_comp);         % Correlation coefficient
        pSm = zeros(num_comp, size(Se, 2));

        for component_ind = 1:num_comp
              
            % Compute spatial correlation with reference
            p = corrcoef(Se(component_ind, :), reference);
            cc(component_ind) = p(1, 2);
            if cc(component_ind)>0
                 pS1 = zscore(Se(component_ind, :));
                 pS1(pS1<1.5)=0; % denoising
                 pSm(component_ind,:)=pS1;   
            else
                 pS1 = -zscore(Se(component_ind, :));
                 pS1(pS1<1.5)=0;
                 pSm(component_ind,:)=pS1;  
            end 
        end
        [c_max, id_select] = max(abs(cc));
        SM_interest(r, m,:) = pSm(id_select, :);
        c_com(r,m)=c_max;
    end
end
[~,id_run]=max(c_com,[],1);
for m = 1:size(SM, 2)
    SM_best(m,:)=SM_interest(id_run(m), m,:);
end
end

