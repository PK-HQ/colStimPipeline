function [mdl]=plotRobustFit(varargin)
%USAGE
%plotRobustFit(x,y,sigTest=[1|0],LineColor='b',LineStyle='--',MarkerStyle='v')

x=varargin{1};
y=varargin{2};

if length(varargin)>=3
    plotSignificance=varargin{3};
else
    plotSignificance=1;
end

if length(varargin)>=4
    color=varargin{4};
else
    color=[1 1 1];
end

if length(varargin)>=5
    linestyle=varargin{5};
else
    linestyle='-';
end

if length(varargin)>=6
    marker=varargin{6};
else
    marker='o';
end

dotSize=50;
dotEdgeColor=[0 0 0];
lineWidth=1.5;

% fit line to data
mdl=fitlm(x,y,'RobustOpts','on');
coeffs=mdl.Coefficients.Estimate;
betaSignificance=mdl.Coefficients.pValue(2) < 0.05;

% generate fitted points for line plotting
xFit=min(x):(max(x)-min(x))/10:max(x);
yFit=coeffs(1)+coeffs(2)*xFit;

if betaSignificance==1 && plotSignificance==1
    lineColor=[0 0 0];
elseif betaSignificance==0 && plotSignificance==1
    lineColor=[.5 .5 .5];
elseif plotSignificance==0
    lineColor=color;
end

h1=scatter(x,y,dotSize,color,'filled','MarkerEdgeColor',dotEdgeColor,'Marker',marker,'LineWidth',1); hold on
h2=plot(xFit,yFit,'Color',lineColor,'LineWidth',lineWidth,'LineStyle',linestyle,'LineWidth',2); hold on

if sum(isnan(xFit))>0
    nanIdx=find(isnan(xFit));
    xFit(nanIdx)=[];
    yFit(nanIdx)=[];
end
    
%{
P=polyfit(x,y,2);
xfit=nanmin(x):abs(nanmin(x):nanmax(x))/100:nanmax(x);
yfit=polyval(P,xfit);
hold on
plot(xfit,yfit,'-','Color',[0.5 0.5 0.5])
%}
%{
bInitial = [0;0;0];
model=@(b,x) b(1) + b(2).*x + b(3).*(x.^2);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
mdl = fitnlm(x,y,model,bInitial);%,'Options',opts);
pvalParam=mdl.Coefficients.pValue(2:3);
params = mdl.Coefficients.Estimate;


if any(pvalParam)<0.05
    lineColor=[0 0 0];
else
    lineColor=[0.5 0.5 0.5];
end
plot(x, model(params,x), 'Color',lineColor,'LineWidth',lineWidth);%lineWidth
%}
end
