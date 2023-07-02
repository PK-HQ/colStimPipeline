function [columnarBitmap,VERpca]=getColumnarBitmapV4(dataStruct,gridSize,gammaCorr, sensitivity, desiredOrts, fastSwitch, plotOrNot)
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

%% Get filenames
filenameStruct=generateFilenames(dataStruct);
pdfFilename=filenameStruct.colmapPDF;
mapsFilename=filenameStruct.colmapMat;
%% Load raw and PCA-ed activity for all orientations, and SNR mask
load(filenameStruct.FFTAmp,'DataCond');
load(filenameStruct.Orientation,'Ort','RespCondPCA','Mask','MapOrt','MapAmpOrt','ColorMap'); %contains mask, RespCondPCA
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
        nCols=7;
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
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            greyNaN(VERpca(:,:,ortNo));
            map=VERpca(:,:,ortNo);
            mapSD=nanstd(map,[],'all');
            mapMean=nanmean(map(:));
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            titleL2=['\sigma=' sprintf('%.1d',mapSD), ' \mu' sprintf('%.1d',mapMean)];
            title({'PCA-ed',titleL1},'FontWeight','normal');
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([-2 2]*10^-3)
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        [ax,h]=suplabel('Columnar bitmap generation','t',[.08 .08 .84 .80]);
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
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            if plotCounter==1
                title({'Normalized',titleL1},'FontWeight','normal');
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
                title({'Equalization',[sprintf('%0g',Ort(ortNo)) '\circ']},'FontWeight','normal');
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
        
        %% 5 Adaptive (local) threshold of the gamma map to obtain bitmap
        plotIdx=[1, ((nRows*nCols)/2 +1)]+4;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            % do adaptive thresholding
            threshold = adaptthresh(double(histEqCol(:,:,ortNo)), sensitivity,...
                'NeighborhoodSize',2*floor(size(histEqCol(:,:,ortNo))/16)+1);
            columnarBitmap(:,:,ortNo) = double(imbinarize(double(histEqCol(:,:,ortNo)),threshold)).*ROImask; %zero mask

            %calculate pixel density
            bitmapNaN=columnarBitmap(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(bitmapNaN(:))*100/nansum(ROImaskNaN(:));

            %plot
            greyNaN(columnarBitmap(:,:,ortNo).*ROImask.*ROItuningMap,.7);
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            titleL2=['PD=' sprintf('%.1f%%',pixelDensity)];
            if plotCounter==1
                title({'Adapt. threshold',titleL1,titleL2},'FontWeight','normal');
            else
                title({titleL1,titleL2},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1]);
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)

        
        
        %% 4 Gamma correction
        plotIdx=[1, ((nRows*nCols)/2 +1)]+3;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            gammCorrImg=imadjust(columnarBitmap(:,:,ortNo),[],[],gammaCorr);
            gammaCol(:,:,ortNo)=gammCorrImg;
            greyNaN(gammaCol(:,:,ortNo).*ROItuningMap);
            
            
            gammaColNaN=gammaCol(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(gammaColNaN(:))*100/nansum(ROImaskNaN(:));
            
            titleL2=['PD=' sprintf('%.1f%%',abs(pixelDensity))];
            if plotCounter==1
                title({'Gamma',[sprintf('%0g',Ort(ortNo)) '\circ'],titleL2},'FontWeight','normal');
            else
                title({[sprintf('%0g',Ort(ortNo)) '\circ'],titleL2},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis('auto')   
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)
        
        
        
        %% 6 Difference bitmaps
        plotIdx=[1, ((nRows*nCols)/2 +1)]+5;plotCounter=1;
        for ortNo=selectedOrt(plotCounter)
            axes(hAx(plotIdx(plotCounter)));
            % do adaptive thresholding
            diffBitmap(:,:,ortNo)=double(gammaCol(:,:,1)>gammaCol(:,:,7));
            %calculate pixel density
            diffbitmapNaN=diffBitmap(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(diffbitmapNaN(:))*100/nansum(ROImaskNaN(:));
            %plot
            greyNaN(diffBitmap(:,:,ortNo).*ROImask.*ROItuningMap,.7);
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            titleL2=['PD=' sprintf('%.1f%%',abs(pixelDensity))];
            if plotCounter==1
                title({'Diff. bitmap',titleL1,titleL2},'FontWeight','normal');
            else
                title({titleL1,titleL2},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1]);
            plotCounter=plotCounter+1;
        end
        
        for ortNo=selectedOrt(plotCounter)
            axes(hAx(plotIdx(plotCounter)));
            % do adaptive thresholding
            diffBitmap(:,:,ortNo)=double(gammaCol(:,:,7)>gammaCol(:,:,1));
            %calculate pixel density
            diffbitmapNaN=diffBitmap(:,:,ortNo).*ROImaskNaN;
            pixelDensity=nansum(diffbitmapNaN(:))*100/nansum(ROImaskNaN(:));
            %plot
            greyNaN(diffBitmap(:,:,ortNo).*ROItuningMap,.7);
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            titleL2=['PD=' sprintf('%.1f%%',abs(pixelDensity))];
            if plotCounter==1
                title({'Diff. bitmap',titleL1,titleL2},'FontWeight','normal');
            else
                title({titleL1,titleL2},'FontWeight','normal');
            end
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis([0 1]);
            plotCounter=plotCounter+1;
        end
        upFontSize(14,.015)
        set(h,'FontSize',16)

        
        %% 7 Overlay with columnar tuning map
        plotIdx=[1, ((nRows*nCols)/2 +1)]+6;plotCounter=1;
        for ortNo=selectedOrt
            axes(hAx(plotIdx(plotCounter)));
            %plot
            greyNaN(diffBitmap(:,:,ortNo).*ROItuningMap,.7);
            titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
            if plotCounter==1
                title({'Sanity',titleL1},'FontWeight','normal');
            else
                title({titleL1},'FontWeight','normal');
            end
            
            axis square;
            addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
            colormap(gray);colorbar;caxis('auto');
            plotCounter=plotCounter+1;
        end
        axes('Position',[.96,0.065,0.01,.85]);
        imagesc(permute(flipdim(ColorMap,1),[1,3,2]));
        set(gca,'YAxisLocation','Right', ...
                'XTick',[], ...
                'YTick',[1,9.5,18], ...
                'YTickLabel',{'180','90','0'});
        ylabel('Orientation ^o', ...
               'FontSize',10, ...
               'FontWeight','Bold');
        upFontSize(14,.015)
        set(h,'FontSize',16)
        export_fig(pdfFilename,'-pdf','-nocrop');

        diffBitmap=diffBitmap(:,:,selectedOrt);
        
        %% Save maps
        save(mapsFilename,'diffBitmap','columnarBitmap','gammaCol','histEqCol','VERpca','normVERpca');

        if plotOrNot==0
            set(groot,'defaultFigureVisible','on')
        end
end
end