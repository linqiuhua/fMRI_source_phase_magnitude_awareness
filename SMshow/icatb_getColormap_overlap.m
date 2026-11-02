function cm = icatb_getColormap(imageValues,useStructural, displayType)
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
            cm(129:256,:) = [[1*ones(20,1);0*ones(108,1)],[(215/255)*ones(20,1);(127/255)*ones(108,1)],0*ones(128,1)];            
            cm(1:128,:) = [(175/255)*ones(128,1),(35/255)*ones(128,1),(222/255)*ones(128,1)];
           
        else
            cm = yellowV; %blue;
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