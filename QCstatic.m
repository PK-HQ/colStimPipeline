function [blurEstimate,gfpStaticMax,mcherryStaticMax]=QCstatic(dataStruct, plotOrNot)

%% What it does
% Loads file specified in dataStruct.field 
% Run whatever you want
%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Get filenames
filenameStruct=generateFilenames(dataStruct); %generates filenames for monkey/run/session in dataStruct
pdfFilename=filenameStruct.neurometricPDF; % here is the save filename

%% Load selected variables from filename
%init
[gfpStatic,mcherryStatic,greenStatic]=deal(nan(512,512,1));

%load if exists
if isfile(filenameStruct.gfpStatic)
    gfpStatic=imread(filenameStruct.gfpStatic);
end
if isfile(filenameStruct.mcherryStatic)
    mcherryStatic=imread(filenameStruct.mcherryStatic);
end
if isfile(filenameStruct.greenStatic)
    greenStatic=imread(filenameStruct.greenStatic);
end

blurEstimate=99; %placeholder
gfpStaticMax=max(double(gfpStatic(:))*100./255);
mcherryStaticMax=max(double(mcherryStatic(:))*100./255);
%sprintf('[%s] GFP: %.1f%% | mCherry: %.1f%%', dataStruct.date, gfpStaticMax, mcherryStaticMax)

figure
[hA,~]=tight_subplot(1,3); 
imgX=[1 size(greenStatic,1)];
imgY=[1 size(greenStatic,2)];

axes(hA(1))
imagesc(double(greenStatic)./255);colorbar;colormap(gray); axis square
title('Green image','FontWeight','Normal')
addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),1,1,3);
caxis([0 1])

axes(hA(2))
imagesc(double(gfpStatic)./255);colorbar;colormap(gray); axis square
title(['GFP static (Peak: ' num2str(gfpStaticMax,'%0.1f%%') ')'],'FontWeight','Normal')
addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),2,1,3);
caxis([0 0.06])

axes(hA(3))
imagesc(double(mcherryStatic)./255);colorbar;colormap(gray); axis square
title(['mCherry static (Peak: ' num2str(mcherryStaticMax,'%0.1f%%') ')'],'FontWeight','Normal')
addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),2,1,3);
caxis([0 .6])
upFontSize(14,.015)

[~,h]=suplabel('Quality of static images','t',[.08 .08 .84 .71]);
set(h,'FontSize',16)

export_fig(pdfFilename,'-pdf','-nocrop','-append');
end