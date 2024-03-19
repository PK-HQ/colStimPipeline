function cax(varargin)
if length(varargin)==0
    modifier=0.9;
    caxisLimits=caxis;caxisLimitsNew=max(abs(caxisLimits));
elseif length(varargin)==1
    modifier=varargin{1};
    caxisLimits=caxis;caxisLimitsNew=max(abs(caxisLimits));
elseif length(varargin)==2
    modifier=varargin{1};
    caxisLimits=[min(varargin{2},[],'all') max(varargin{2},[],'all')];
    caxisLimitsNew=max(abs([min(varargin{2},[],'all') max(varargin{2},[],'all')]));
end
colorbar
caxis('auto')


if  round(caxisLimits(1),2)<0 & caxisLimits(2)>0 % + - range
    caxis([-caxisLimitsNew caxisLimitsNew]*modifier);
    colormap(fireice)
elseif round(caxisLimits(1),2)>=0 & caxisLimits(2)>0 % + range
    caxis([0 caxisLimitsNew]*modifier); 
    colormap(fire)
elseif caxisLimits(1)<0 & round(caxisLimits(2),2)<=0 % - range
    caxis([-caxisLimitsNew 0]*modifier);
    colormap(ice)
end
end