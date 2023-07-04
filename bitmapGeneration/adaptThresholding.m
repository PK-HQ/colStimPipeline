function [diffBitmap,pixelDensity,estDutycycle]=adaptThresholding(gammaCol,Sensitivity,ROImask,ROImaskNaN)

% Get bitmap for each set of columns
for ortNo=1:size(gammaCol,3)
      % do adaptive thresholding
    threshold = adaptthresh(double(gammaCol(:,:,ortNo)),Sensitivity,...
        'NeighborhoodSize',2*floor(size(gammaCol(:,:,ortNo))/16)+1);
    columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold)).*ROImask; %zero mask
    
    
    % calculate duty cycle and pixel density
    bitmapNaN=columnarBitmap(:,:,ortNo).*ROImaskNaN;
    pixelDensity(ortNo)=nansum(bitmapNaN(:))*100/nansum(ROImaskNaN(:));
    estDutycycle(ortNo)=sqrt(pixelDensity*2/100)*5/10;
    
end

% Subtract them to eliminate overlap, these are the final bitmaps
diffBitmap(:,:,1)=double(columnarBitmap(:,:,1)>columnarBitmap(:,:,7));
diffBitmap(:,:,2)=double(columnarBitmap(:,:,7)>columnarBitmap(:,:,1));


end