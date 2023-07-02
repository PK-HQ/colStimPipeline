function [diffBitmap,VERpca,PCAExplTotal,pixelDensities,bitmapEnergies]=getColumnarBitmapV4(dataStructRef, dataStructSession,gridSize,gammaCorr, sensitivity, desiredOrts, fastSwitch, plotOrNot, isSummary)
%% Change log
% 1. Suspect that there are pixels with abnormally high value in 90-165deg patterns, try solve by clipping outliers 
% outside +-3sd at [-1 ,1]

% 2. F0 footprint crop by using defining ROI per orientation map that is > 70% amplitude, then crop all to that ROI

%% What it does
% Loads raw and PCA-ed responses and mask. 
% Collects SNR and ROI masked blank subtracted VER and PCA-ed VER (VER,
% VERpca)
% 1. SNR and ROI masking 
% 2. Rescale + histogram equalization, i.e. histEqCol
% 3. Gamma correction, i.e. gammaCol
% 4. Thresholding, i.e. columnarBitmap

%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Define column selection method
columnSelectionMethod='MapOrt';
orangeLightPD=56.1;

%% Get filenames
if isSummary
    filenameStructSession=generateFilenames(dataStructSession);
    load(filenameStructSession.TS)
end

filenameStructRef=generateFilenames(dataStructRef);
pdfFilename=filenameStructSession.neurometricPDF;
mapsFilename=filenameStructSession.colmapMat;

%% Load raw and PCA-ed activity for all orientations, and SNR mask
load(filenameStructRef.FFTAmp,'DataCond');
load(filenameStructRef.Orientation,'Ort','RespCondPCA','Mask','MapOrt','MapAmpOrt','ColorMap','nPCAComp','PCAExpl'); %contains mask, RespCondPCA
%conds
nOrt=numel(Ort); %excl blanks
if fastSwitch
    selectedOrt=find(ismember(Ort,desiredOrts)==1);
else
    selectedOrt=1:nOrt;
end
%ROI
imgDims=size(RespCondPCA(:,:,1));
imgX=[0 imgDims(1)];
imgY=[0 imgDims(2)];
ROIx=[1 imgDims(1)];%[round(.016/.016) imgDims(1)];
ROIy=[1 imgDims(2)];%[round(3.2/.016) imgDims(2)];
ROIMask=NaN(imgX(2),imgY(2));
ROIMask(ROIy(1):ROIy(2),ROIx(1):ROIx(2))=1;

%n
nBlanks=size(DataCond,3)-nOrt;
DataCondBlankSubt=DataCond(:,:,nBlanks+1:nOrt+nBlanks)-mean(DataCond(:,:,1:nBlanks),3);

%% Extract masked VER and pca-ed VER
%ROI mask
ROImask=double(Mask);
ROImaskNaN=ROImask;
ROImaskNaN(ROImaskNaN==0)=NaN;
ROItuningMap=ROImask.*MapAmpOrt;

%% Figure params
switch fastSwitch %only processes 0 and 90 columns (for speed)
    case 1
        nRows=2;
        nCols=6;
        for ortNo=selectedOrt
            %VER for orientation, masked with SNR mask
            VERraw(:,:,ortNo)=DataCondBlankSubt(:,:,ortNo).*ROImaskNaN;
            VERpca(:,:,ortNo)=RespCondPCA(:,:,ortNo).*ROImaskNaN;

            %VER for orientation, masked with ROI mask (to isolate site)
            VERraw(:,:,ortNo)=VERraw(:,:,ortNo).*ROIMask;
            VERpca(:,:,ortNo)=VERpca(:,:,ortNo).*ROIMask;
        end

        %% 1 Visualize PCA-ed data
        figure('name','PCA-ed responses')
        colormap(gray);
        [hAx,~]=tight_subplot(nRows,nCols,[.025 .025]);
        plotIdx=[1, ((nRows*nCols)/2 +1)];plotCounter=1;
        
        PCAExplTotal=sum(PCAExpl(1:nPCAComp));
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            greyNaN(VERpca(:,:,ortNo));
            map=VERpca(:,:,ortNo);
            mapSD=nanstd(map,[],'all');
            mapMean=nanmean(map(:));
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            if ortNo==1
                title({['PCA-ed (nComp=' num2str(nPCAComp,'%.0f') ')'],['Expl. var=' num2str(PCAExplTotal,'%.0f') '%'],titleL1},'FontWeight','normal');
            else
                xlabel('X (mm)');
                ylabel('Y (mm)');

            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([-2 2]*10^-3)
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        [~,h]=suplabel('Columnar bitmap generation','t',[.08 .08 .84 .80]);
        set(h,'FontSize',16)
        
        %% 2 Normalize it
        colormap(gray);
        plotIdx=[1, ((nRows*nCols)/2 +1)]+1;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            %rescale to norm
            map=VERpca(:,:,ortNo);
            normMap=rescale(map,0,1); 
            %clips image, remap to -1 to 1
            normMap = normMap .* ROImaskNaN;
            normVERpca(:,:,ortNo)=normMap;

            %plot
            greyNaN(normVERpca(:,:,ortNo).*ROItuningMap);
            if plotCounter==1
                title({'Normalized'},'FontWeight','normal');
            else
                title({titleL1},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;
            plotCounter=plotCounter+1;
        end

        %% 3 Histogram equalization to maximize dynamic range across image
        plotIdx=[1, ((nRows*nCols)/2 +1)]+2;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            histEqCol(:,:,ortNo)=adapthisteq(normVERpca(:,:,ortNo),'NumTiles',imgDims./gridSize,'Range','full');
            greyNaN(histEqCol(:,:,ortNo).*ROItuningMap);colorbar;colormap(fireice)
            if plotCounter==1
                title({'Equalization'},'FontWeight','normal');
            else
                title({[sprintf('%0g',Ort(ortNo)) '\circ']},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis('auto')
            plotCounter=plotCounter+1;
            offwarning
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)

        %% 4 Gamma correction
        plotIdx=[1, ((nRows*nCols)/2 +1)]+3;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            gammCorrImg=imadjust(histEqCol(:,:,ortNo),[],[],gammaCorr);
            gammaCol(:,:,ortNo)=gammCorrImg;
            greyNaN(gammaCol(:,:,ortNo).*ROItuningMap);
            if plotCounter==1
                title({'Gamma corr.'},'FontWeight','normal');
            else
                title({[sprintf('%0g',Ort(ortNo)) '\circ']},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis('auto')   
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)

        %% 5 Adaptive (local) threshold of the gamma map to obtain bitmap
        plotIdx=[1, ((nRows*nCols)/2 +1)]+4;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            % do adaptive thresholding
            threshold = adaptthresh(double(gammaCol(:,:,ortNo)), sensitivity,...
                'NeighborhoodSize',2*floor(size(gammaCol(:,:,ortNo))/16)+1);
            columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold)).*ROImask; %zero mask

            %calculate pixel density
            bitmapNaN=columnarBitmap(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(bitmapNaN(:))*100/nansum(ROImaskNaN(:));
            equivalentDutyCycle=sqrt(pixelDensity*2/100)*5/10;
            %plot
            greyNaN(columnarBitmap(:,:,ortNo).*ROImask.*ROItuningMap,.7);
            titleL1=['DC=' sprintf('%.2f',equivalentDutyCycle) ' (' sprintf('%.1f%%',pixelDensity) ')'];
            if plotCounter==1
                title({'Adapt. threshold',titleL1},'FontWeight','normal');
            else
                title({titleL1},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1]);
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)

        %% 6 Difference bitmaps
        plotIdx=[1, ((nRows*nCols)/2 +1)]+5;plotCounter=1;
        for ortNo=selectedOrt(plotCounter)
            axes(hAx(plotIdx(plotCounter)));
            % do adaptive thresholding
            diffBitmap(:,:,ortNo)=double(columnarBitmap(:,:,1)>columnarBitmap(:,:,7));
            %calculate pixel density
            diffbitmapNaN=diffBitmap(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(diffbitmapNaN(:))*100/nansum(ROImaskNaN(:));
            equivalentDutyCycle=sqrt(pixelDensity*2/100)*5/10;
            
            if isSummary
              %calculate bitmap energy
              bitmapEnergy=calculateBitmapEnergy(TS, orangeLightPD, pixelDensity);
              pixelDensities(ortNo)=pixelDensity;
              bitmapEnergies(ortNo)=bitmapEnergy;
            else
              pixelDensities(ortNo)=pixelDensity;
              bitmapEnergies(ortNo)=NaN;
            end

            %plot
            greyNaN(diffBitmap(:,:,ortNo).*ROImask.*ROItuningMap,.7);
            titleL1=['DC=' sprintf('%.2f',equivalentDutyCycle) ' (' sprintf('%.1f%%',pixelDensity) ')'];
            if plotCounter==1
                title({'Diff. bitmap',titleL1},'FontWeight','normal');
            else
                title({titleL1},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1]);
            plotCounter=plotCounter+1;
        end
        
        for ortNo=selectedOrt(plotCounter)
            axes(hAx(plotIdx(plotCounter)));
            % do adaptive thresholding
            diffBitmap(:,:,ortNo)=double(columnarBitmap(:,:,7)>columnarBitmap(:,:,1));
            %calculate pixel density
            diffbitmapNaN=diffBitmap(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(diffbitmapNaN(:))*100/nansum(ROImaskNaN(:));
            equivalentDutyCycle=sqrt(pixelDensity*2/100)*5/10;
            
            if isSummary
              %calculate bitmap energy
              bitmapEnergy=calculateBitmapEnergy(TS, orangeLightPD, pixelDensity);
              pixelDensities(ortNo)=pixelDensity;
              bitmapEnergies(ortNo)=bitmapEnergy;
            else
              pixelDensities(ortNo)=pixelDensity;
              bitmapEnergies(ortNo)=NaN;
            end
            
            %plot
            greyNaN(diffBitmap(:,:,ortNo).*ROItuningMap,.7);
            titleL1=['DC=' sprintf('%.2f',equivalentDutyCycle) ' (' sprintf('%.1f%%',pixelDensity) ')'];
            if plotCounter==1
                title({'Diff. bitmap',titleL1},'FontWeight','normal');
            else
                title({titleL1},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1]);
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)

        export_fig(pdfFilename,'-pdf','-nocrop');
       
        diffBitmap(:,:,selectedOrt);
        bitmapEnergies=bitmapEnergies(selectedOrt);
        
        
        
        %% Save maps
        save(mapsFilename,'diffBitmap','columnarBitmap','gammaCol','histEqCol','VERpca','normVERpca');

        %% Show plot and goldilocks plot (for iterating sensitivity)
        %{
        figure('name','Threshold')

        [~,OrtSubset]=find(Ort==[0 90]');

        for ortNo=1:OrtSubset
            subplot(2,2,ortNo);

            imgDist=gammaCol(:,:,ortNo);

            % do adaptive thresholding
            sensitivity=0.01;
            threshold = adaptthresh(double(gammaCol(:,:,ortNo)), sensitivity);
            columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold));

            BWnan=columnarBitmap(:,:,ortNo).*ROImaskNaN;
            greyNaN(columnarBitmap(:,:,ortNo).*ROImask,.7);

            pixelDensity=nansum(BWnan(:))/numel(ROImaskNaN);

            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            titleL2=['PD=' sprintf('%.1d',pixelDensity)];
            title({titleL1,titleL2},'FontWeight','normal');
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1])%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);
        end
        suplabel(['Visual-evoked response (rescaled-PCA + HE + gamma-corrected + threshold, Threshold=' num2str(threshPrctile) '%)'],'t');
        export_fig(pdfFilename,'-pdf','-nocrop','-append');
        %}

        if plotOrNot==0
            set(groot,'defaultFigureVisible','on')
        end
end
end