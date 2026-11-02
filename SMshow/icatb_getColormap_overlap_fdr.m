function cm = icatb_getColormap(imageValues,useStructural, displayType)
%--------------------------------------------------------------------------
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
%       imageValues = 1 pos neg
%                     2 pos
%       useStructural = 1 if you are overlaying on a structural image
%
%   OUTPUT:
%       cm = colormap
%
% #############################################################

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
            cm(129:256,:) = [[1*ones(4,1); 0.9*ones(31,1); 0.8*ones(91,1); 0*ones(2,1)], ...
                 [(215/255)*ones(4,1); 0.4*ones(31,1); 0.1*ones(91,1); (142/255)*ones(2,1)], ...
                 [0*ones(4,1); 0.8*ones(31,1); 0.1*ones(91,1); 0*ones(2,1)]];
            cm(1:128,:) = [[(175/255)*ones(112,1);(1)*ones(16,1)],...
                [(35/255)*ones(112,1);(0.5)*ones(16,1)],...
                [(222/255)*ones(112,1);(0.2)*ones(16,1)]];
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