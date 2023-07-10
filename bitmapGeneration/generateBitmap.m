function [VERpca,normVERpca,histEqCol,gammaCol,columnarBitmap,columnarBitmapStats,optostimBitmap]=generateBitmap(selectedOrt,responsePCAImages,ROImask,bitmapParams)
%% README
% Function processes PCA-ed columnar responses and outputs optostim bitmaps
% for biasing experiments that apply optostim as manipulation. Steps are:
% 1. Restrict processing to desired area via ROI mask
% 2. Enhance columnar responses by normalizing to [0, 1], then equalizing
%    contrast across entire image (grid-wise) with histogram equalization
% 3. Isolate column centers with gamma correction
% 4. Select columns with adaptive thresholding, and subtract to ensure
% bitmaps are non-overlapping
% 
% 
% 



%% Make ROI mask
ROImaskNaN=createNanMask(ROImask);

% get wante
responseImages=responsePCAImages(:,:,selectedOrt);
imgDims=size(responseImages(:,:,1));

for ortNo=1:size(responseImages,3)
    
    %% Mask pca response with SNR mask (selected in RunDA, typically d'>=8 or RMS>=0.01)
    VERpca(:,:,ortNo)=responseImages(:,:,ortNo).*ROImaskNaN;
    VERpca(:,:,ortNo)=VERpca(:,:,ortNo);

    %% normalize to [0,1]
    normVERpca(:,:,ortNo)=rescale(VERpca(:,:,ortNo),0,1).* ROImaskNaN; 
  
    %% histogram equalization,
    histEqCol(:,:,ortNo)=adapthisteq(normVERpca(:,:,ortNo),'NumTiles',imgDims./bitmapParams.gridSize,'Range','full');

    
    %% gamma corr
    gammaCol(:,:,ortNo)=imadjust(histEqCol(:,:,ortNo),[],[],bitmapParams.gammaCorrFactor);
    
    
     %% do adaptive thresholding
    nbhdSize=2*floor(size(gammaCol(:,:,ortNo))/16)+1; % 2*(64/16) +1
    threshold = adaptthresh(double(gammaCol(:,:,ortNo)),bitmapParams.sensitivity,...
        'NeighborhoodSize',nbhdSize);
    columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold)).*ROImask; %zero mask

    
    %% calculate bitmap statistocs: duty cycle and pixel density
    bitmapNaN=columnarBitmap(:,:,ortNo).*ROImaskNaN;
    columnarBitmapStats.pixelDensity(ortNo)=nansum(bitmapNaN(:))*100/nansum(ROImaskNaN(:));
    columnarBitmapStats.estDutycycle(ortNo)=sqrt(columnarBitmapStats.pixelDensity(ortNo)*2/100)*5/10;
  
    
    
end

% Subtract them to eliminate overlap, these are the final bitmaps
optostimBitmap(:,:,1)=double(columnarBitmap(:,:,1)>columnarBitmap(:,:,2));
optostimBitmap(:,:,2)=double(columnarBitmap(:,:,2)>columnarBitmap(:,:,1));

end