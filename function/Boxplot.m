function Boxplot(data_sz, data_hc, pv)
% BoxPlot - Plot boxplots for a given feature in schizophrenia (SZ) and healthy control (HC) groups.
%
% Inputs:
%   - data_sz: feature vector of SZ group [K1 x 1] 
%            K1 is the number of schizophrenia patients.
%   - data_hc: feature vector of HC group [K2 x 1] 
%            K2 is the number of healthy controls.
%   - pv: p-value of a two sample t-test between SZ and HC groups
%
%% construct group labels
groupLabels = [repmat({'SZ'}, numel(data_sz), 1); repmat({'HC'}, numel(data_hc), 1)];

%% Set boxplot parameters and custom styles
h = boxplot([data_sz;data_hc], groupLabels, 'Widths', 0.6, 'Colors', 'k', 'Symbol', '+');
hold on;
box_colors = [141/255,206/255,250/255; 1, 165/255, 195/255];
box_h = findobj(gca, 'Tag', 'Box');
for j = 1:length(box_h)
    patch(get(box_h(j), 'XData'), get(box_h(j), 'YData'), box_colors(j, :), ...
        'FaceAlpha', 0.6, 'EdgeColor', 'none'); % Ìí¼Ó°ëÍ¸Ã÷Ìî³ä
end

box_h = findobj(gca, 'Tag', 'Box');
median_h = findobj(gca, 'Tag', 'Median');
outlier_h = findobj(gca, 'Tag', 'Outliers');

set(box_h(2), 'Color', [1, 1, 1]);
set(median_h(2), 'Color', [197,22,22]/255, 'LineWidth', 1);
set(outlier_h(1), 'MarkerEdgeColor', [197,22,22]/255);

set(box_h(1), 'Color', [1,1,1]);
set(median_h(1), 'Color', [0,100,160]/255, 'LineWidth', 1);
set(outlier_h(2), 'MarkerEdgeColor', [0,100,160]/255);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
ylim([min(min(data_sz),min(data_hc))-(max(max(data_sz),max(data_hc))-min(min(data_sz),max(max(data_sz),max(data_hc))))*0.2 max(max(data_sz),max(data_hc))*1.5])
box on

%% Annotate statistical significance if p < 0.05
if pv < 0.05
yl = ylim;
y_star = yl(2) - 0.15 * range(yl);
text(1.5, y_star, '*', 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
