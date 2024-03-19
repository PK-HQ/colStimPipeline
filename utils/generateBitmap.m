function [VERpca,normVERpca,histEqCol,gammaCol,columnarBitmap,bitmapData]=generateBitmap(selectedOrt,imagingData,bitmapData,blockID)
%% README
% Function processes PCA-ed columnar responses and outputs optostim bitmaps
% for biasing experiments that apply optostim as manipulation. Steps are:
% 1. Restrict processing to desired area via ROI mask
% 2. Enhance columnar responses by normalizing to [0, 1], then equalizing
%    contrast across entire image (grid-wise) with histogram equalization
% 3. Isolate column centers with gamma correction
% 4. Select columns with adaptive thresholding, and subtract to ensure
% bitmaps are non-overlapping


%% Make ROI mask
ROImaskNaN=createNanMask(imagingData.mask(:,:,blockID));

% get wanted
responseImages=imagingData.ortpca(:,:,selectedOrt,blockID);
imgDims=size(responseImages(:,:,1));

for ortNo=1:size(responseImages,3)
    
    %% Mask pca response with SNR mask (selected in RunDA, typically d'>=8 or RMS>=0.01)
    VERpca(:,:,ortNo)=responseImages(:,:,ortNo).*ROImaskNaN;
    %VERpca(:,:,ortNo)=VERpca(:,:,ortNo);

    %% normalize to [0,1]
    normVERpca(:,:,ortNo)=rescale(VERpca(:,:,ortNo),0,1).* ROImaskNaN; 
  
    %% histogram equalization,
    histEqCol(:,:,ortNo)=adapthisteq(normVERpca(:,:,ortNo),'NumTiles',imgDims./bitmapData.gridSize(:,blockID),'Range','full');
    offwarning;
    
    %% gamma corr
    gammaCol(:,:,ortNo)=imadjust(histEqCol(:,:,ortNo),[],[],bitmapData.gammaCorrFactor(1,ortNo,blockID));
     
     %% do adaptive thresholding
    nbhdSize=2*floor(size(gammaCol(:,:,ortNo))/16)+1; % 2*(64/16) +1
    threshold = adaptthresh(double(gammaCol(:,:,ortNo)),bitmapData.sensitivity(:,blockID),...
        'NeighborhoodSize',nbhdSize);
    columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold)).*imagingData.mask(:,:,blockID); %zero mask

    %% calculate bitmap statistocs: duty cycle and pixel density
    bitmapNaN=columnarBitmap(:,:,ortNo).*ROImaskNaN;
    bitmapData.pixDensity(ortNo,blockID)=nansum(bitmapNaN(:))*100/nansum(ROImaskNaN(:));
    bitmapData.DC(ortNo,blockID)=sqrt(bitmapData.pixDensity(ortNo,blockID)*2/100)*5/10;
      
end

% Subtract them to eliminate overlap, these are the final bitmaps
optostimBitmap(:,:,1)=double(columnarBitmap(:,:,1)>columnarBitmap(:,:,2));
optostimBitmap(:,:,2)=double(columnarBitmap(:,:,2)>columnarBitmap(:,:,1));
bitmapData.columnarbitmap(:,:,:,blockID)=optostimBitmap;
end