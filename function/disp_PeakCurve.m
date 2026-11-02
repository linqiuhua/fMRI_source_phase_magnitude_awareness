function disp_PeakCurve(submap_mag,submap_phase,location_x,location_y,location_z,com_id,sr)
% disp_PeakCurve - Display spatial magnitude and phase curves around 
%                region peak voxel in x,y,z directions
%
% disp_PeakCurve - Display magnitude and phase polarity curves along three spatial directions 
%                  centered at a peak voxel.
% Inputs:
%   - submap_mag: magnitude map within sub-region [X x Y x Z]
%            X, Y, Z represent three spatial dimensions.
%   - submap_phase: phase map within sub-region [X x Y x Z]
%   - location_x/y/z: coordinate of peak voxel
%   - com_id: component ID for plot formatting
%   - sr: sub-region ID for plot formatting
%
% Output:
%   Generate 2x3 figure panel showing:
%   - top row: magnitude curves with phase transition markers
%   - bottom row: phase polarity curves

%% x-direction
figure;subplot(2,3,1);
magnitude = squeeze(submap_mag(:,location_y,location_z));
phase = squeeze(submap_phase(:, location_y, location_z));
[~, peak_idx] = max(magnitude);
phase_sign = sign(phase);phase_sign(phase == 0) = 0;

% detect phase polarity transitions
polar_changes=[];
for ii=1:length(phase_sign)-1
    if phase_sign(ii)*phase_sign(ii+1)<0
        polar_changes=[polar_changes ii];
    end
end
if numel(polar_changes)~=0
polar_changes(polar_changes == max(find(phase_sign))) = [];
end
% annotate phase transitions with directional markers
plot(magnitude, 'Color', [1, 0.647, 0], 'LineWidth', 3);hold on;
for i = 1:length(polar_changes)
    change_idx = polar_changes(i);

    % negative-to-positive transition (red upward arrow)
    % positive-to-negative transition (blue downward arrow)
    if change_idx < peak_idx
        if phase_sign(change_idx) > 0
            plot(change_idx, magnitude(change_idx), 'rv', 'MarkerFaceColor',  [0, 112/255, 192/255], 'MarkerEdgeColor',  [0, 112/255, 192/255], 'MarkerSize', 8);            
        else
            plot(change_idx, magnitude(change_idx), 'g^', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor',  [1, 0, 0], 'MarkerSize', 8);
        end
    elseif change_idx >= peak_idx
        if phase_sign(change_idx) > 0            
            plot(change_idx, magnitude(change_idx), 'g^', 'MarkerFaceColor',   [1, 0, 0], 'MarkerEdgeColor',   [1, 0, 0], 'MarkerSize', 8);
        else
            plot(change_idx, magnitude(change_idx), 'rv', 'MarkerFaceColor',  [0, 112/255, 192/255], 'MarkerEdgeColor', [0, 112/255, 192/255], 'MarkerSize', 8);
        end
    end
end


    % adaptively ajust
    if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end;
        set(gca,'FontName','Arial','FontSize',16);
        title(['Location: [',num2str(location_x),',',num2str(location_y),',',num2str(location_z),']']);
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end; 
 % end
set(gca,'FontName','Arial','FontSize',16);
title('x-direction');

% display phase polarity
subplot(2,3,4);
tmp = sign(squeeze(submap_phase(:, location_y, location_z)));tmp(tmp < 0) = 0;
stem(tmp, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 2, 'Color', 'r');
hold on;tmp = sign(squeeze(submap_phase(:, location_y, location_z)));tmp(tmp > 0) = 0;
stem(tmp, 'Marker', 'o', 'MarkerFaceColor', [0,112/255,192/255], 'MarkerEdgeColor', [0,112/255,192/255], ...
'LineStyle', '-', 'LineWidth', 2, 'Color', [0,112/255,192/255]);
tmp = sign(squeeze(submap_phase(:, location_y, location_z)));tmp(tmp ~= 0) = NaN;
stem(find(~isnan(tmp)), tmp(~isnan(tmp)), 'Marker', 'o', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'LineStyle', 'none', 'LineWidth', 0.5);
[maxv,maxidx]=(max(squeeze(submap_mag(:,location_y,location_z))));
hold on;plot(maxidx,sign(squeeze(submap_phase(maxidx,location_y,location_z))),'o', 'MarkerFaceColor', [1, 0.647, 0], 'Color', [1, 0.647, 0]);
set(gca,'FontName','Arial','FontSize',16);

% set plot parameters
% if com_id==1 && sr==1 
%     xlim([10 45]); ylim([-2 2]);set(gca, 'YTick', []);
% elseif com_id==2 && sr==1 
%     xlim([0 37]); ylim([-2 2]);set(gca, 'YTick', []);
% elseif com_id==3 && sr==3 
%     xlim([20 53]); ylim([-2 2]);set(gca, 'YTick', []);
% else

% adaptively ajust
if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
    xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end;
    set(gca,'FontName','Arial','FontSize',16);
    if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
    xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end; 
    ylim([-2 2]);set(gca, 'YTick', []);
% end


%% y-direction
subplot(2,3,2);
magnitude = squeeze(submap_mag(location_x,:,location_z));
phase = squeeze(submap_phase(location_x,:,location_z));
[~, peak_idx] = max(magnitude);
phase_sign = sign(phase);phase_sign(phase == 0) = 0;
polar_changes=[];
for ii=1:length(phase_sign)-1
    if phase_sign(ii)*phase_sign(ii+1)<0
        polar_changes=[polar_changes ii];
    end
end

if numel(polar_changes)~=0
polar_changes(polar_changes == max(find(phase_sign))) = [];
end
subplot(2,3,2);plot(magnitude, 'Color', [1, 0.647, 0], 'LineWidth', 3);hold on;
for i = 1:length(polar_changes)
    change_idx = polar_changes(i);
    if change_idx < peak_idx
        if phase_sign(change_idx) > 0
            plot(change_idx, magnitude(change_idx), 'rv', 'MarkerFaceColor',  [0, 112/255, 192/255], 'MarkerEdgeColor',  [0, 112/255, 192/255], 'MarkerSize', 8);
            
        else
            plot(change_idx, magnitude(change_idx),  'g^', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor',  [1, 0, 0], 'MarkerSize', 8);
        end
    elseif change_idx >= peak_idx
        if phase_sign(change_idx) > 0
            plot(change_idx, magnitude(change_idx), 'g^', 'MarkerFaceColor',  [1, 0, 0], 'MarkerEdgeColor', [1, 0, 0], 'MarkerSize', 8);
        else
            plot(change_idx, magnitude(change_idx), 'rv', 'MarkerFaceColor',  [0, 112/255, 192/255], 'MarkerEdgeColor',  [0, 112/255, 192/255], 'MarkerSize', 8);
        end
    end
end
 
% if com_id==1 && sr==1 
%     xlim([41 62]); ylim([0 5]);
%     elseif com_id==2 && sr==1 
%         xlim([14 51]); ylim([0 10]);
%     elseif com_id==3 && sr==3 
%         xlim([30 63]); ylim([0 7]);
%     else

    if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end;
        set(gca,'FontName','Arial','FontSize',16);
        title(['Location: [',num2str(location_x),',',num2str(location_y),',',num2str(location_z),']']);
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end; 
% end
set(gca,'FontName','Arial','FontSize',16);
title('y-direction');


subplot(2,3,5);tmp = sign(squeeze(submap_phase(location_x,:,location_z)));tmp(tmp < 0) = 0;
stem(tmp, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 2, 'Color', 'r');
hold on;tmp = sign(squeeze(submap_phase(location_x,:,location_z)));
tmp(tmp > 0) = 0;    stem(tmp, 'Marker', 'o', 'MarkerFaceColor', [0,112/255,192/255], 'MarkerEdgeColor', [0,112/255,192/255], ...
'LineStyle', '-', 'LineWidth', 2, 'Color', [0,112/255,192/255]);
tmp = sign(squeeze(submap_phase(location_x,:,location_z)));tmp(tmp ~= 0) = NaN;
stem(find(~isnan(tmp)), tmp(~isnan(tmp)), 'Marker', 'o', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'LineStyle', 'none', 'LineWidth', 0.5);

[maxv,maxidx]=(max(squeeze(submap_mag(location_x,:,location_z))));
hold on;plot(maxidx,sign(squeeze(submap_phase(location_x,maxidx,location_z))),'o', 'MarkerFaceColor', [1, 0.647, 0], 'Color', [1, 0.647, 0]);   
set(gca,'FontName','Arial','FontSize',16);
% if com_id==1 && sr==1 
%     xlim([41 62]); ylim([-2 2]);set(gca, 'YTick', []);
%     elseif com_id==2 && sr==1 
%         xlim([14 51]); ylim([-2 2]);set(gca, 'YTick', []);
%     elseif com_id==3 && sr==3 
%         xlim([30 63]); ylim([-2 2]);set(gca, 'YTick', []);
%     else
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end;
        set(gca,'FontName','Arial','FontSize',16);
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end; 
        ylim([-2 2]);set(gca, 'YTick', []);
% end


%% z-direction
subplot(2,3,3);    
magnitude = squeeze(submap_mag(location_x,location_y,:));
phase = squeeze(submap_phase(location_x,location_y,:));
[peak_value, peak_idx] = max(magnitude);
phase_sign = sign(phase);phase_sign(phase == 0) = 0;

polar_changes=[];
for ii=1:length(phase_sign)-1
    if phase_sign(ii)*phase_sign(ii+1)<0
        polar_changes=[polar_changes ii];
    end
end
if numel(polar_changes)~=0
polar_changes(polar_changes == max(find(phase_sign))) = [];
end
subplot(2,3,3);plot(magnitude, 'Color', [1, 0.647, 0], 'LineWidth', 3);hold on;
for i = 1:length(polar_changes)
    change_idx = polar_changes(i);
    if change_idx < peak_idx
        if phase_sign(change_idx) > 0
            plot(change_idx, magnitude(change_idx), 'g^', 'MarkerFaceColor',  [1, 0, 0], 'MarkerEdgeColor',  [1, 0, 0], 'MarkerSize', 10);            
        else
            plot(change_idx, magnitude(change_idx),  'rv', 'MarkerFaceColor',  [0, 112/255, 192/255], 'MarkerEdgeColor',  [0, 112/255, 192/255], 'MarkerSize', 8);
        end
    elseif change_idx >= peak_idx
        if phase_sign(change_idx) > 0
            plot(change_idx, magnitude(change_idx),  'g^', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor',  [1, 0, 0], 'MarkerSize', 8);
        else
            plot(change_idx, magnitude(change_idx), 'rv', 'MarkerFaceColor',  [0, 112/255, 192/255], 'MarkerEdgeColor',  [0, 112/255, 192/255], 'MarkerSize', 8);
        end
    end
end

% if com_id==1 && sr==1 
%     xlim([5 35]); ylim([0 5]);
%     elseif com_id==2 && sr==1 
%         xlim([8 45]); ylim([0 10]);
%     elseif com_id==3 && sr==3 
%         xlim([7 45]); ylim([0 7]);
%     else

    if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end;
        set(gca,'FontName','Arial','FontSize',16);
        title(['Location: [',num2str(location_x),',',num2str(location_y),',',num2str(location_z),']']);
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end; 
% end
set(gca,'FontName','Arial','FontSize',16);
title('z-direction');


subplot(2,3,6); tmp = sign(squeeze(submap_phase(location_x,location_y,:)));tmp(tmp < 0) = 0;
stem(tmp, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 2, 'Color', 'r');
hold on;tmp = sign(squeeze(submap_phase(location_x,location_y,:)));
tmp(tmp > 0) = 0;
stem(tmp, 'Marker', 'o', 'MarkerFaceColor', [0,112/255,192/255], 'MarkerEdgeColor', [0,112/255,192/255], ...
'LineStyle', '-', 'LineWidth', 2, 'Color', [0,112/255,192/255]);    
tmp = sign(squeeze(submap_phase(location_x,location_y,:)));tmp(tmp ~= 0) = NaN;
stem(find(~isnan(tmp)), tmp(~isnan(tmp)), 'Marker', 'o', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'LineStyle', 'none', 'LineWidth', 0.5);
if length(find(squeeze(submap_mag(location_x,location_y,:))))~=0
xlim([min(find(squeeze(submap_mag(location_x,location_y,:))))-5 max(find(squeeze(submap_mag(location_x,location_y,:))))+5]);end;ylim([-2 2]);set(gca, 'YTick', []);
[maxv,maxidx]=(max(squeeze(submap_mag(location_x,location_y,:))));
hold on;plot(maxidx,sign(squeeze(submap_phase(location_x,location_y,maxidx))),'o', 'MarkerFaceColor', [1, 0.647, 0], 'Color', [1, 0.647, 0]);
set(gca,'FontName','Arial','FontSize',16);

% if com_id==1 && sr==1 
%     xlim([5 35]); ylim([-2 2]);set(gca, 'YTick', []);
%     elseif com_id==2 && sr==1 
%         xlim([8 45]); ylim([-2 2]);set(gca, 'YTick', []);
%     elseif com_id==3 && sr==3 
%         xlim([7 45]); ylim([-2 2]);set(gca, 'YTick', []);
%     else
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end;
        set(gca,'FontName','Arial','FontSize',16);
        if length(find(squeeze(submap_mag(:,location_y,location_z))))~=0
        xlim([min(find(squeeze(submap_mag(:,location_y,location_z))))-5 max(find(squeeze(submap_mag(:,location_y,location_z))))+5]);end; 
        ylim([-2 2]);set(gca, 'YTick', []);
% end

end