function [b0, b1, b2, a_b, n_b, e] = piecewise_linear_fit(am, ntxyz, fig_op,sub_ind,showinfo)
% piecewise_linear_fit - Adaptive piecewise linear fitting
%
% Inputs:
%   - am: mean activation magnitude [1 x N]
%            N = number of magnitude voxel cubes
%   - ntxyz: mean number of polarity transitions along the x,y,z directions [1 x N]
%   - fig_op: flag to enable (1) or disable (0) result plotting.
%   - showinfo: figure title string

% Outputs:
%   - b0: slope of first segment.
%   - b1: intercept of first segment.
%   - b2: slope of second segment.
%   - a_b: optimal breakpoint of am.
%   - n_b: ntxyz value at optimal breakpoint.
%   - e: minimum mean square error (MSE) of optimal fit.

%% Breakpoint optimization
SEARCH_PRECISION = 0.001;  % Step size for breakpoint search

if max(ntxyz) == 0
    b0 = 0; b1 = 0; b2 = 0;     
    a_b = 0;  n_b = 0; e=0;
    return;
end

ii = 1;
for a_b_candidate = min(am):SEARCH_PRECISION:(max(am) - SEARCH_PRECISION)

    % ---- split data ----
    idx1 = (am <= a_b_candidate);
    idx2 = (am >  a_b_candidate);
    if nnz(idx1) < 2 || nnz(idx2) < 2
        e_search(ii) = Inf;
        ii = ii + 1;
        continue;
    end

    % fit first segment
    warning('off');
    p1 = polyfit(am(idx1), ntxyz(idx1), 1);

    % Evaluate fitting results at the candidate breakpoint
    if (max(am(idx1)) == 0 || isnan(p1(1)))
        n_b = 0;
    else
        n_b = polyval(p1, a_b_candidate);
    end

    % fit second segment with fixed starting point (connect at breakpoint)
    [fit_segment1, ~] = constrained_linear_fit(am(idx2) - a_b_candidate, ntxyz(idx2), n_b);
    b2 = fit_segment1.p1;

    % compute combined MSE
    fit_segment2 = n_b + b2 * (am(idx2) - a_b_candidate);

    e_search(ii) = (sum((ntxyz(idx1) - polyval(p1, am(idx1))).^2) + ...
                    sum((ntxyz(idx2) - fit_segment2).^2)) / ...
                   (nnz(idx1) + nnz(idx2));

    ii = ii + 1;
end


% optimal fit determination
[e, ii_optimal] = min(e_search); 
a_b = min(am) + (SEARCH_PRECISION * (ii_optimal - 1)); % optimal breakpoint

%% final parameter estimation
fit_segment1_optimal = polyfit(am(am <= a_b), ntxyz(am <= a_b), 1);
n_b = polyval(fit_segment1_optimal, a_b); 
b1=fit_segment1_optimal(1);b0=fit_segment1_optimal(2);
[fit_segment2_optimal, ~] = constrained_linear_fit(am(am > a_b)- a_b, ntxyz(am > a_b), n_b);
b2 = fit_segment2_optimal.p1; 

%% result visualization
if fig_op
    figure;
    % ---- scatter ----
    x1 = am(am > a_b);  y1 = ntxyz(am > a_b);
    x2 = am(am <= a_b); y2 = ntxyz(am <= a_b);
    if sub_ind <= 40
        mkFace = [141,206,250]/255;
        lnCol  = [4, 99, 176]/255;
    else
        mkFace = [255,165,195]/255;
        lnCol  = [200, 0, 0]/255;
    end
    plot(x1, y1, 'o', 'MarkerFaceColor', mkFace, ...
         'MarkerEdgeColor', 'none', 'MarkerSize', 5);
    hold on;
    plot(x2, y2, 'o', 'MarkerFaceColor', mkFace, ...
         'MarkerEdgeColor', 'none', 'MarkerSize', 5);

    % ---- fitted lines ----
    xFit1 = linspace(min(am), a_b, 200);
    yFit1 = polyval(fit_segment1_optimal, xFit1);
    xFit2 = linspace(a_b, max(am), 200);
    yFit2 = n_b + b2 * (xFit2 - a_b);
    plot(xFit1, yFit1, '--', 'Color', lnCol, 'Linewidth', 1.5);
    plot(xFit2, yFit2, '--', 'Color', lnCol, 'Linewidth', 1.5);
    plot(a_b, n_b, 'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', lnCol, ...
    'MarkerEdgeColor', lnCol);

    % ---- axis style ----
    set(gca, 'FontName', 'Arial', 'FontSize', 10);
    axis equal;
    box on;
    xlim([min(x2)-0.4 max(x1)+0.4]);
    ylim([min(y2)-0.1 max(y1)+0.2]);
    xlabel('a_m', 'Interpreter','tex', 'FontSize', 14);
    ylabel('n_{txyz}', 'Interpreter','tex', 'FontSize', 14);

    infoStr = sprintf(['$b_0$=%.2f,  $b_1$=%.2f,  $b_2$=%.2f\n' ...
                       '$a_b$=%.2f,  $n_b$=%.2f,  $e$=%.2f'], ...
                       b0, b1, b2, a_b, n_b, e);

    text(0.5, 0.97, infoStr, ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'Interpreter','latex', ...
        'FontName','Arial', ...
        'FontSize',16, ...
        'Color','k');

    titleStr = sprintf('PWL for %s', showinfo);
    title(titleStr, 'Interpreter','none', 'FontName','Arial', 'FontSize', 10);


    hold off;
end

end

function [fit_result, gof] = constrained_linear_fit(x, y, fixed_intercept)
% Linear regression with fixed intercept

% Inputs:
%   - x, y: data for fitting
%   - fixed_intercept: constrained y-intercept value

% Outputs:
%   - fit_result: fit object containing slope parameter
%   - gof: A structure containing goodness of fit info

    [xData, yData] = prepareCurveData(x, y);

    % Define the fit model
    func = [num2str(fixed_intercept), '+ p1 * x'];
    ft = fittype(func, 'independent', 'x', 'dependent', 'y');

    % Set up fit options
    opts = fitoptions(ft);
    opts.StartPoint = [0.5]; 
    % opts.Lower = 0; 
    % opts.Upper = Inf; 

    % Perform the fit
    [fit_result, gof] = fit(xData, yData, ft, opts);
end
% 
% 
% function [b0, b1, b2, a_b, n_b, e] = piecewise_linear_fit(am, ntxyz, fig_op,sub_ind)
% % piecewise_linear_fit - Adaptive piecewise linear fitting
% %
% % Inputs:
% %   - am: mean activation magnitude [1 x N]
% %            N = number of magnitude voxel cubes
% %   - ntxyz: mean number of polarity transitions along the x,y,z directions [1 x N]
% %   - fig_op: flag to enable (1) or disable (0) result plotting.
% %
% % Outputs:
% %   - b0: slope of first segment.
% %   - b1: intercept of first segment.
% %   - b2: slope of second segment.
% %   - a_b: optimal breakpoint of am.
% %   - n_b: ntxyz value at optimal breakpoint.
% %   - e: minimum mean square error (MSE) of optimal fit.
% 
% %% Breakpoint optimization
% SEARCH_PRECISION = 0.001;  % Step size for breakpoint search
% 
% if max(ntxyz) == 0
%     b0 = 0; b1 = 0; b2 = 0;     
%     a_b = 0;  n_b = 0;    
%     return;
% end
% 
% ii = 1;
% for a_b_candidate = min(am):SEARCH_PRECISION:(max(am) - SEARCH_PRECISION)        
%     % fit first segment
%     warning('off'); 
%     p1 = polyfit(am(am <= a_b_candidate), ntxyz(am <= a_b_candidate), 1);
% 
%     % Evaluate fitting results at the candidate breakpoint
%     if (max(am(am <= a_b_candidate)) == 0 || isnan(p1(1)))
%         n_b = 0; 
%     else
%         n_b = polyval(p1, a_b_candidate); 
%     end
% 
%     % fit second segment with fixed starting point (connect at breakpoint)
%     [fit_segment1, ~] = constrained_linear_fit(am(am > a_b_candidate)...
%                 - a_b_candidate, ntxyz(am > a_b_candidate), n_b); 
%     b2 = fit_segment1.p1;
% 
%     % compute combined MSE
%     fit_segment2 = n_b + fit_segment1.p1 * (am(am > a_b_candidate) - a_b_candidate);
%     e_search(ii) = (sum((ntxyz(am <= a_b_candidate) - polyval(p1, am(am <= a_b_candidate))).^2)...
%                     + sum((ntxyz(am > a_b_candidate) - fit_segment2).^2)) / ...
%                    (length(am(am <= a_b_candidate)) + length(am(am > a_b_candidate)));
%     ii = ii + 1;
% end
% 
% % optimal fit determination
% [e, ii_optimal] = min(e_search); 
% a_b = min(am) + (SEARCH_PRECISION * (ii_optimal - 1)); % optimal breakpoint
% 
% %% final parameter estimation
% fit_segment1_optimal = polyfit(am(am <= a_b), ntxyz(am <= a_b), 1);
% n_b = polyval(fit_segment1_optimal, a_b); 
% b1=fit_segment1_optimal(1);b0=fit_segment1_optimal(2);
% [fit_segment2_optimal, ~] = constrained_linear_fit(am(am > a_b)- a_b, ntxyz(am > a_b), n_b);
% b2 = fit_segment2_optimal.p1; 
% 
% %% result visualization
% if fig_op
%  if sub_ind <=40
%     figure;
%     x1 = am(am > a_b);y1 = ntxyz(am > a_b);
%     plot(x1, y1, 'o', 'MarkerFaceColor', [141,206,250]/255, ...
%          'MarkerEdgeColor', 'none', 'MarkerSize', 5);
%     hold on;
%     x2 = am(am <= a_b);y2 = ntxyz(am <= a_b);
%     plot(x2, y2, 'o', 'MarkerFaceColor', [141,206,250]/255, ...
%          'MarkerEdgeColor', 'none', 'MarkerSize', 5);    
%     plot([0:0.01:a_b], polyval(fit_segment1_optimal, [0:0.01:a_b]), ...
%          'Color', [4, 99, 176]/255, 'Linewidth', 1.5, 'LineStyle', '--');
%     plot(a_b:0.01:max(am), n_b + b2*((a_b:0.01:max(am))-a_b), ...
%          'Color', [4, 99, 176]/255, 'Linewidth', 1.5, 'LineStyle', '--');
%     set(gca, 'FontName', 'Arial', 'FontSize', 12);
%     axis equal;
%      hold off;
%  else
%     figure;
%     x1 = am(am > a_b);y1 = ntxyz(am > a_b);
%     plot(x1, y1, 'o', 'MarkerFaceColor', [255,165,195]/255, ...
%          'MarkerEdgeColor', 'none', 'MarkerSize', 5);
%     hold on;
%     x2 = am(am <= a_b);y2 = ntxyz(am <= a_b);
%     plot(x2, y2, 'o', 'MarkerFaceColor', [255,165,195]/255, ...
%          'MarkerEdgeColor', 'none', 'MarkerSize', 5);    
%     plot([0:0.01:a_b], polyval(fit_segment1_optimal, [0:0.01:a_b]), ...
%          'Color', [200, 0, 0]/255, 'Linewidth', 1.5,  'LineStyle', '--');
%     plot(a_b:0.01:max(am), n_b + b2*((a_b:0.01:max(am))-a_b), ...
%          'Color', [200, 0, 0]/255, 'Linewidth', 1.5,  'LineStyle', '--');
%     set(gca, 'FontName', 'Arial', 'FontSize', 12);   
%     axis equal;
% 
%     hold off;
%   end
% end
% end
% 
% function [fit_result, gof] = constrained_linear_fit(x, y, fixed_intercept)
% % Linear regression with fixed intercept
% 
% % Inputs:
% %   - x, y: data for fitting
% %   - fixed_intercept: constrained y-intercept value
% 
% % Outputs:
% %   - fit_result: fit object containing slope parameter
% %   - gof: A structure containing goodness of fit info
% 
%     [xData, yData] = prepareCurveData(x, y);
% 
%     % Define the fit model
%     func = [num2str(fixed_intercept), '+ p1 * x'];
%     ft = fittype(func, 'independent', 'x', 'dependent', 'y');
% 
%     % Set up fit options
%     opts = fitoptions(ft);
%     opts.StartPoint = [0.5]; 
% 
%     % Perform the fit
%     [fit_result, gof] = fit(xData, yData, ft, opts);
% end
