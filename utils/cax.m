function cax(varargin)
% varargin{1:2} = modifier, img
if length(varargin)==0
    modifier=0.8;
    caxisLimits=caxis;caxisLimitsNew=max(abs(caxisLimits));
elseif length(varargin)==1
    modifier=varargin{1};
    caxisLimits=caxis;caxisLimitsNew=max(abs(caxisLimits));
elseif length(varargin)==2
    modifier=varargin{1};
    caxisLimits=[min(varargin{2},[],'all') max(varargin{2},[],'all')];
    caxisLimitsNew=max(abs([min(varargin{2},[],'all') max(varargin{2},[],'all')]));
end
colorbarHandle=colorbar;
colorbarHandle.FontSize = 12; % Set the desired font size
caxis('auto')
if  round(caxisLimits(1),2,'significant')<0 & caxisLimits(2)>0 % + - range
    caxis([-caxisLimitsNew caxisLimitsNew]*modifier);
    cmap = fireice;             % Get the current colormap
    cmap(end+1, :) = [1 1 1]; % Append gray to colormap
    %colormap(cmap);              % Set the modified colormap
    % Set the color of NaNs to use the last color in colormap
    set(gca, 'Color', cmap(end, :));
elseif round(caxisLimits(1),2,'significant')>=0 & caxisLimits(2)>0 % + range
    caxis([0 caxisLimitsNew]*modifier); 
    %colormap(fire)
elseif caxisLimits(1)<0 & round(caxisLimits(2),2,'significant')<=0 % - range
    caxis([-caxisLimitsNew 0]*modifier);
    %colormap(ice)
end
end