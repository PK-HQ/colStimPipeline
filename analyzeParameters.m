%% Load data

% def mainPath
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
load([mainPath 'Chip/Meta/summary/neurometrics20230413.mat'])
pdfFilename=[mainPath 'Chip/Meta/summary/summary20230413.mat'];

% Outlier removal
outlierRemoval=0;

%% Define biasing and predictors
% Delta biasing (V, H, V-H)
usableSessionIdx=find(abs(deltaSortedCont(:,2)-0.5)<0.5); %sessions with BL delta < 0.5 +- 0.2
vDeltas=deltaSortedCont(usableSessionIdx,3)-deltaSortedCont(usableSessionIdx,2); %deltaSortedCont(:,2)
hDeltas=deltaSortedCont(usableSessionIdx,1)-deltaSortedCont(usableSessionIdx,2);
vhDeltas=deltaSortedCont(usableSessionIdx,3)-deltaSortedCont(usableSessionIdx,1); %deltaSortedCont(:,2)

% Define predictors (V, H, V-H)
predictors_x_V=[gfpStaticMaxCont(usableSessionIdx)';
    mcherryStaticMaxCont(usableSessionIdx)';
    PCAExplTotalCont(usableSessionIdx)';
    pixelDensitiesCont(usableSessionIdx,2)';
    nBlobsCont(usableSessionIdx,2)';
    medianBlobAreasCont(usableSessionIdx,2)';
    powerDensitiesCont(usableSessionIdx,2)';
    visualrCont(usableSessionIdx,2)';
    visualrRatiosCont(usableSessionIdx,2)';
    PCArFullCont(usableSessionIdx,4)';
    PCArFullCont(usableSessionIdx,3)';
    PCArRatiosCont(usableSessionIdx,2)';
    optoIntensityCont(usableSessionIdx,2)'];

predictors_xLabel_V={'GFP intensity (%)';
    'mCherry intensity (%)';
    'PCA expl. var. (%)';
    'Pixel density (%)';
    'No. of columns';
    'Median blob area';
    'Bitmap energy (a.u.)';
    'Correlation visual';
    'Correlation ratio visual';
    'Correlation PCA (Con.)';
    'Correlation PCA (Incon.)';
    'Correlation ratio PCA (Incon/Con)';
    'Intensity visual'};


predictors_x_H=[gfpStaticMaxCont(usableSessionIdx)';
    mcherryStaticMaxCont(usableSessionIdx)';
    PCAExplTotalCont(usableSessionIdx)';
    pixelDensitiesCont(usableSessionIdx,1)';
    nBlobsCont(usableSessionIdx,1)';
    medianBlobAreasCont(usableSessionIdx,1)';
    powerDensitiesCont(usableSessionIdx,1)';
    visualrCont(usableSessionIdx,1)';
    visualrRatiosCont(usableSessionIdx,1)';
    PCArFullCont(usableSessionIdx,1)';
    PCArFullCont(usableSessionIdx,2)';
    PCArRatiosCont(usableSessionIdx,1)';
    optoIntensityCont(usableSessionIdx,1)'];

predictors_xLabel_H={'GFP intensity (%)';
    'mCherry intensity (%)';
    'PCA expl. var. (%)';
    'Pixel density (%)';
    'No. of columns';
    'Median blob area';
    'Bitmap energy (a.u.)';
    'Correlation visual';
    'Correlation ratio visual';
    'Correlation PCA (Con.)';
    'Correlation PCA (Incon.)';
    'Correlation ratio PCA (Incon/Con)';
    'Intensity visual'};

predictors_x_VH=[gfpStaticMaxCont(usableSessionIdx)';
    mcherryStaticMaxCont(usableSessionIdx)';
    PCAExplTotalCont(usableSessionIdx)';
    pixelDensitiesCont(usableSessionIdx,2)'-pixelDensitiesCont(usableSessionIdx,1)';
    nBlobsCont(usableSessionIdx,2)'-nBlobsCont(usableSessionIdx,1)';
    medianBlobAreasCont(usableSessionIdx,2)'-medianBlobAreasCont(usableSessionIdx,1)';
    powerDensitiesCont(usableSessionIdx,2)'-powerDensitiesCont(usableSessionIdx,1)';
    visualrCont(usableSessionIdx,2)'-visualrCont(usableSessionIdx,1)';
    visualrRatiosCont(usableSessionIdx,2)'-visualrRatiosCont(usableSessionIdx,1)';
    PCArFullCont(usableSessionIdx,3)'-PCArFullCont(usableSessionIdx,2)';
    PCArFullCont(usableSessionIdx,4)'-PCArFullCont(usableSessionIdx,1)';
    PCArRatiosCont(usableSessionIdx,2)'-PCArRatiosCont(usableSessionIdx,1)';
    optoIntensityCont(usableSessionIdx,2)'-optoIntensityCont(usableSessionIdx,1)'];

predictors_xLabel_VH={'\Delta GFP intensity (%)';
    '\Delta mCherry intensity (%)';
    '\Delta PCA expl. var. (%)';
    '\Delta pixel density (%)';
    '\Delta no. of columns';
    '\Delta median blob area';
    '\Delta bitmap energy (a.u.)';
    '\Delta Correlation visual';
    '\Delta Correlation ratio visual';
    '\Delta Correlation PCA (Diff Con.)';
    '\Delta Correlation PCA (Diff Incon.)';
    '\Delta Correlation ratio PCA (Diff Incon/Con)';
    '\Delta Intensity visual'};

%% Linear regression (Bias and parameters, V-H)
figure('name','Session-wise PCA corr.')
[hAx,~]=tight_subplot(2,1,[.1 .25],[],[.2 .1]);%delete(hAx([4 12]));

axes(hAx(1))
visColors={'b','r'}; 
optoColors={'b','r'};
optoMarkers={'v','^'};
linestyles={'--','--'};
lineWidth=2;
dotSize=50;
for optoCond=[1 2]
    % def X and Y
    x=1:size(deltaSortedCont,1);
    
    % just to adjust for the indexing in deltaSortedCont
    if optoCond==1
        optoCondIdx=1;
    elseif optoCond==2
        optoCondIdx=3;
    end
    
    y=deltaSortedCont(usableSessionIdx,optoCondIdx)-deltaSortedCont(usableSessionIdx,2);
    h2=plot(x,y,'Color',visColors{optoCond},'LineWidth',lineWidth,'LineStyle',linestyles{optoCond},'LineWidth',2); hold on
    h1=scatter(x,y,dotSize,optoColors{optoCond},'MarkerEdgeColor',optoColors{optoCond},'Marker',optoMarkers{optoCond},'LineWidth',1); hold on
    
    % option for outlier removal
    switch outlierRemoval
        case 1
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
    end

    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end
yLabel='\beta';
xLabel='Session';
ylabel(yLabel);
xlabel(xLabel);
% labels
addSkippedTicks(0,16,1,'x')
ylim([-.75 .75]);addSkippedTicks(-.75,.75,.25,'y');
title('Psychometric: \beta');
yline(0,'k--')
upFontSize(14,0.005);


axes(hAx(2))
visColors={'b', 'r', 'b', 'r'}; 
optoColors={'b', 'b', 'r', 'r'};
optoMarkers={'v','v','^','^'};
linestyles={'--',':',':','--'};
lineWidth=2;
dotSize=50;
for optoCond=[1 4]%1:size(PCArFullCont,2)
    % def X and Y
    x=1:size(PCArFullCont,1);
    y=PCArFullCont(:,optoCond);
    h2=plot(x,y,'Color',visColors{optoCond},'LineWidth',lineWidth,'LineStyle',linestyles{optoCond},'LineWidth',2); hold on
    h1=scatter(x,y,dotSize,optoColors{optoCond},'MarkerEdgeColor',optoColors{optoCond},'Marker',optoMarkers{optoCond},'LineWidth',1); hold on
    
    % option for outlier removal
    switch outlierRemoval
        case 1
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
    end

    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end

title('Neurometric:{\it r}_{PCA}','Interpreter','tex') %Low contrast versus PCA reference
yLabel='Correlation';
xLabel='Session';
ylabel(yLabel);
xlabel(xLabel);
% labels
addSkippedTicks(0,16,1,'x')
ylim([-.75 .75]);addSkippedTicks(-.75,.75,.25,'y');
yline(0,'k--')
upFontSize(14,0.005);
%suplabel('Session-wise psychometric and neurometric data','t',[.08 .08 .84 .78]);
export_fig(pdfFilename,'-pdf','-nocrop');





figure('name','Session-wise PCA corr.')
[hAx,~]=tight_subplot(2,1,[.2 .25]);%delete(hAx([4 12]))


axes(hAx(1))
visColors={[138,43,226]*.9/255}; 
optoColors={[138,43,226]*.9/255};
optoMarkers={'diamond'};
linestyles={'--'};
lineWidth=2;
dotSize=50;
for optoCond=1
    % def X and Y
    x=1:size(deltaSortedCont,1);
        
    y=deltaSortedCont(usableSessionIdx,3)-deltaSortedCont(usableSessionIdx,1);
    h2=plot(x,y,'Color',visColors{optoCond},'LineWidth',lineWidth,'LineStyle',linestyles{optoCond},'LineWidth',2); hold on
    h1=scatter(x,y,dotSize,optoColors{optoCond},'MarkerEdgeColor',optoColors{optoCond},'Marker',optoMarkers{optoCond},'LineWidth',1); hold on
    
    % option for outlier removal
    switch outlierRemoval
        case 1
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
    end

    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end
title('Psychometric: \Delta \beta');
xLabel='Session';
yLabel='\Delta \beta';
ylabel(yLabel);
xlabel(xLabel);
% labels
addSkippedTicks(0,16,1,'x')
ylim([-.75 .75]);addSkippedTicks(-.75,.75,.25,'y');
yline(0,'k--')
upFontSize(14,0.005);


axes(hAx(2))
visColors={[138,43,226]*.9/255,[138,43,226]*.9/255}; 
optoColors={[138,43,226]*.9/255,[138,43,226]*.9/255};
optoMarkers={'diamond','diamond'};
linestyles={'--',':'};
lineWidth=2;
dotSize=50;
for optoCond=[1 2]%1:size(PCArFullCont,2)
    % def X and Y
    x=1:size(PCArFullCont,1);
    

        % just to adjust for the indexing in deltaSortedCont
    if optoCond==1
        y=PCArFullCont(:,4)-PCArFullCont(:,1);
    elseif optoCond==2
        optoCondIdx=PCArFullCont(:,3)-PCArFullCont(:,2);
    end
    
    h2=plot(x,y,'Color',visColors{optoCond},'LineWidth',lineWidth,'LineStyle',linestyles{optoCond},'LineWidth',2); hold on
    h1=scatter(x,y,dotSize,optoColors{optoCond},'MarkerEdgeColor',optoColors{optoCond},'Marker',optoMarkers{optoCond},'LineWidth',1); hold on
    
    % option for outlier removal
    switch outlierRemoval
        case 1
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
    end

    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end
yline(0,'k--')
title('Neurometric: \Delta Low contrast versus PCA reference')
yLabel='\Delta Correlation';
xLabel='Session';
ylabel(yLabel);
xlabel(xLabel);
% labels
addSkippedTicks(0,16,1,'x')
ylim([-1 1]);addSkippedTicks(-1,1,.25,'y');
upFontSize(14,0.005);
xlim([0 16])
suplabel('Session-wise \Delta psychometric and neurometric data','t',[.08 .08 .84 .78]);
export_fig(pdfFilename,'-pdf','-nocrop','-append');
















%% Linear regression (Bias and parameters, V-H)
figure('name','Regression, V-optostim')
[hAx,~]=tight_subplot(3,5,[.125 .05]);%delete(hAx([4 12]));
imagingQualityIdx=1:3;
optostimParameters=6:9;
optostimQualityIdx=11:15;
plotIdx=[imagingQualityIdx optostimParameters optostimQualityIdx]; delete(hAx([4:5 10]))

for measure=1:size(predictors_x_V,1)-1
    y=vDeltas;
    x=predictors_x_V(measure,:);
    
    if outlierRemoval==1
        [~,idxX]=rmoutliers(x);
        idx=(idxX')>0;
        x=x(~idx);
        y=y(~idx);
    elseif outlierRemoval==2
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
        x=x(~idx);
        y=y(~idx);
    end
    
    xLabel=predictors_xLabel_V{measure};
    axes(hAx(plotIdx(measure)));
    
    if intersect(plotIdx(measure),imagingQualityIdx)
        c=[34,139,34]/255;
    elseif intersect(plotIdx(measure),optostimParameters)
        c=[255, 191, 0]/255;
    elseif intersect(plotIdx(measure),optostimQualityIdx)
        c=[115, 41, 149]/255;
    end
    
    plotRobustFit(x,y,1,c)
    upFontSize(14,0.005);axis square
    
    % labels
    if plotIdx(measure)==1
        yLabel='V bias (\Delta\beta_{V})';
        ylabel(yLabel);
    end

    yMin=roundData(min(y),0.25);
    yMax=roundData(max(y),0.25);

    %ylim([yMin yMax]);
    addSkippedTicks(yMin,yMax,.25,'y');
    
    xlabel(xLabel);
    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end
1;
suplabel('Read-write parameters on biasing (\Delta\beta_{V})','t',[.08 .08 .84 .76]);
upFontSize(14,0.005);
export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Linear regression (Bias and parameters, H)
figure('name','Regression, H-optostim')
[hAx,~]=tight_subplot(3,5,[.125 .05]);%delete(hAx([4 12]));
imagingQualityIdx=1:3;
optostimParameters=6:9;
optostimQualityIdx=11:15;
plotIdx=[imagingQualityIdx optostimParameters optostimQualityIdx]; delete(hAx([4:5 10]))

for measure=1:size(predictors_x_H,1)-1
    y=hDeltas;
    x=predictors_x_H(measure,:);
    
    if outlierRemoval==1
        [~,idxX]=rmoutliers(x);
        idx=(idxX')>0;
        x=x(~idx);
        y=y(~idx);
    elseif outlierRemoval==2
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
        x=x(~idx);
        y=y(~idx);
    end
    
    xLabel=predictors_xLabel_H{measure};
    axes(hAx(plotIdx(measure)));
    
    if intersect(plotIdx(measure),imagingQualityIdx)
        c=[34,139,34]/255;
    elseif intersect(plotIdx(measure),optostimParameters)
        c=[255, 191, 0]/255;
    elseif intersect(plotIdx(measure),optostimQualityIdx)
        c=[115, 41, 149]/255;
    end
    
    plotRobustFit(x,y,1,c)
    upFontSize(14,0.005);axis square
    
    % labels
    if plotIdx(measure)==1
        yLabel='H bias (\Delta\beta_{H})';
        ylabel(yLabel);
    end

    yMin=roundData(min(y),0.25);
    yMax=roundData(max(y),0.25);
    
    %ylim([yMin yMax]);
    addSkippedTicks(yMin,yMax,.25,'y');
    
    xlabel(xLabel);
    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end
1;
suplabel('Read-write parameters on biasing (\Delta\beta_{H})','t',[.08 .08 .84 .76]);
upFontSize(14,0.005);
export_fig(pdfFilename,'-pdf','-nocrop','-append');


%% Linear regression (Bias and parameters, V-H)
figure('name','Regression, V-H optostim')
[hAx,~]=tight_subplot(3,5,[.125 .05]);%delete(hAx([4 12]));
imagingQualityIdx=1:3;
optostimParameters=6:9;
optostimQualityIdx=11:15;
plotIdx=[imagingQualityIdx optostimParameters optostimQualityIdx]; delete(hAx([4:5 10]))

for measure=1:size(predictors_x_VH,1)-1
    y=vhDeltas;
    x=predictors_x_VH(measure,:);
    
    if outlierRemoval==1
        [~,idxX]=rmoutliers(x);
        idx=(idxX')>0;
        x=x(~idx);
        y=y(~idx);
    elseif outlierRemoval==2
        [~,idxX]=rmoutliers(x);
        [~,idxY]=rmoutliers(y);
        idx=(idxX+idxY')>0;
        x=x(~idx);
        y=y(~idx);
    end
    
    xLabel=predictors_xLabel_VH{measure};
    axes(hAx(plotIdx(measure)));
    
    if intersect(plotIdx(measure),imagingQualityIdx)
        c=[34,139,34]/255;
    elseif intersect(plotIdx(measure),optostimParameters)
        c=[255, 191, 0]/255;
    elseif intersect(plotIdx(measure),optostimQualityIdx)
        c=[115, 41, 149]/255;
    end
    
    plotRobustFit(x,y,1,c)
    upFontSize(14,0.005);axis square
    
    % labels
    if plotIdx(measure)==1
        yLabel='V - H bias (\Delta\beta_{V-H})';
        ylabel(yLabel);
    end

    yMin=roundData(min(y),0.25);
    yMax=roundData(max(y),0.25);

    %ylim([yMin yMax]);
    addSkippedTicks(yMin,yMax,.25,'y');
    
    xlabel(xLabel);
    if any(abs(y)>=100)
        ax=gca; ax.YAxis.Exponent = 3;
    end
end
suplabel('Read-write parameters on biasing (\Delta\beta_{V-H})','t',[.08 .08 .84 .76]);
upFontSize(14,0.005);
export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Multivariate regression
%{
set1=[3, ... %PCA expl
    5 6 7,... %bitmap columns, area, energy
    9 10]; %corr PCA, intensity visual
set2=[1 2 5 6 7 9 10];

mdlVH=stepwiselm(predictors_x_VH(1:end-1,:)',vhDeltas,'linear')

mdlV=stepwiselm(predictors_x_V([1:10],:)',vDeltas,'linear')

mdlH=stepwiselm(predictors_x_H(set2,:)',hDeltas,'linear')


coeffs=mdl.Coefficients.Estimate;
eqnY=coeffs(1)+coeffs(2)*x;
dotSize=25;
dotFillColor=[1 1 1];
dotEdgeColor=[0 0 0];
lineWidth=2;
lineColor=[0 0 0];
%}
figure
scatter(x,eqnY,dotSize,dotFillColor,'filled','MarkerEdgeColor',dotEdgeColor); hold on
plot(x,eqnY,'Color',lineColor,'LineWidth',lineWidth);


%% Heatmap correlation
figure
r = corr([vhDeltas predictors_x_VH'], 'rows','complete');
r=round(tril(r),2);
idx = tril(r);
r(~idx) = nan;
h = heatmap(r, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
labels = ['\Delta bias'; predictors_xLabel_VH];
h.XDisplayLabels = labels;
h.YDisplayLabels = labels; 
colormap(fireice)

