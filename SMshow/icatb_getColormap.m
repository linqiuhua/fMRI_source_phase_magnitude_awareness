function cm = icatb_getColormap(imageValues,useStructural, displayType)
% cm = icatb_getColormap(numOfComp,negValues,useStructural)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: gets colormap based on number of components, whether the data is
%   absoulte value and if there is a structural image. A matlab file specified 
%   by icatb_defaults contains the colorsmaps it uses
%
% #########################################################################
% 
% USAGE:
%
% -cm = icatb_getColormap(numOfComp,negValues,useStructural)
%   INFO: Gets colormap
%   ARGUMENTS:
%       numOfComp = number of components you need a color map for
%       imageValues = 1 pos neg
%                       2 pos
%                       3 abs
%       useStructural = 1 if you are overlaying on a structural image
%
%   OUTPUT:
%       cm = colormap
%
% #############################################################
% 
% LICENSING:
% 
% 
%------------------------------------------------------




%load colorfile defaults
icatb_defaults;
global COLORMAP_FILE;

if ~exist('displayType', 'var')
    displayType = 'other';
end

load(COLORMAP_FILE);

%--set up colormap
colorLength = 64;
  
if(imageValues == 1)
    R=[1 0 0];
      cm(129:256,:) = [1*ones(128,1),[1/128:1/64:0.5,97/192:1/192:1]',0*ones(128,1)];%(1/128:1/128:1)';
    cm(1:128,:) = [[(74/256:-(75-3)/256/128:3/256)';0.01],[(188/256:-(188-79)/256/128:80/256)';0.313],[(238/256:-(238-211)/256/128:212/256)';0.8278;0.827;0.8262;0.8254]];        
else
    cm = yellowV;
end
              
colorbarSkip = size(cm,1)/colorLength;
tempCM = [cm(1:colorbarSkip:end,:)];
cm = tempCM;

if (useStructural)
    gm = linspace(0, 1, 64);
    if size(gm, 1) == 1
        gm = gm';
    end
    structCM = [gm, gm, gm];
    colorbarSkip = size(structCM,1)/colorLength;
    cm=[cm;structCM(1:colorbarSkip:end, :)];
end

colorbarSkip = ceil(size(cm,1)/256);
cm = cm(1:colorbarSkip:end, :);