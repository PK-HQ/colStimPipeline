function [bitmapData]=convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
    currentBlockStruct, conversionType, blockID, ...
    pdfFilename, plotFlag,saveFlagBMP,saveFlag)
disp('Generating bitmap...')

%% Generate bitmap here
for bitmapNo = 1:size(bitmapData.columnarbitmapCoreg,3) % for each input image in camera space
    
    % Define image and contour map wanted
    ort=bitmapData.orts(1,bitmapNo,blockID);
    bitmapCamspace=bitmapData.columnarbitmapCoreg(:,:,bitmapNo,blockID); % grab input image wanted
    contourMaskLevel=bitmapData.gaussianContourLevel(1,bitmapNo,blockID); % grab contour mask desired to extract n-columns
    
    [imagingData, ort, bitmapCamspace, contourMaskLevel]=getColumnSelectionMask(bitmapData, imagingData, blockID, bitmapNo, plotFlag, saveFlag, pdfFilename);

    if ort==0 || ort ==90 % Process only 0 or 90 deg orientations (for speed)
        switch conversionType
            case {'cam2proj'} % Convert from camera space (512 x 512px) to projector space (1920 x 1080px)
               %% Plot 1: Select central n-columns with a contour (from 2D gaussian)
                [gaussianMask, bitmapCamspacePostMask] = isolateColumns(bitmapData, imagingData, blockID, bitmapNo, bitmapCamspace);
                     
               %% Plot 2: Resize to projector space
                [bitmapProjSpaceResized] = resizeToSpace(imagingData, bitmapData, blockID, bitmapCamspacePostMask, conversionType);

               %% Plot 3: Cropping (of y) and adding black bars (of x)
                bitmapProjSpaceBB = addBlackBars(bitmapData, bitmapProjSpaceResized, conversionType);

                %% Apply camera-projector transformation
                bitmapProjSpaceAligned=applyCamProjAlignment(imagingData, blockID, bitmapProjSpaceBB, conversionType);

                %% Save to bitmapData
                % Calculate bitmap stats for columns within the contour (total pixels on, % pixel on, spatial duty cycle)
                nanContour=double(gaussianMask);
                nanContour(nanContour==0)=NaN;
                contourROI=nanContour.*bitmapCamspace;
               [labeledImage,columnAreas,boundaryBox]=getBlobAreas(contourROI);

                % save
                bitmapData.pixelsON(bitmapNo,blockID)=nansum(contourROI(:)>0);
                bitmapData.pixelsONDensity(bitmapNo,blockID)=nansum(contourROI(:)>0)*100/nansum(contourROI(:)==0); %prctPixON=((DC/(5/10))^2)/(2/100);
                bitmapData.columnarbitmapTFprojspace(:,:,bitmapNo,blockID)=bitmapProjSpaceAligned;
                bitmapData.columnarbitmapTFcamspace(:,:,bitmapNo,blockID)=bitmapCamspacePostMask;
                bitmapData.nColumns(bitmapNo,blockID)=numel(columnAreas);
                bitmapData.medianColumnAreas(bitmapNo,blockID)=median(columnAreas);
                
                %% Save columnar optostim BMP
                if ispc & saveFlagBMP==1
                    bmpPath=sprintf('T:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                    bmpFilename = sprintf('O%05gHE%04gG%03gS%05gC%02g.bmp',...
                        ort*100,bitmapData.gridSize,bitmapData.gammaCorrFactor,bitmapData.sensitivity*10000,bitmapData.gaussianContourLevel(bitmapNo));

                    if ~exist(bmpPath, 'dir')
                         mkdir(bmpPath)
                    end

                    % Save the bitmapData
                    disp(bmpFilename)
                    imwrite(bitmapProjSpaceAligned,[bmpPath bmpFilename]);
                end
                
                
                %% Power density calculations
                chamberID='L';
                LEDpercent=30;
                [bitmapData.adjustedSPD_uW(1,bitmapNo,blockID), bitmapData.ledpower_mW(1,bitmapNo,blockID)] = calculateSPD(behavioralData, imagingData, bitmapData,...
                    bitmapNo, chamberID,LEDpercent,blockID,0);

                %% Save the figures
                 [~,h]=suplabel(['Converting bitmap for projector dimensions (' num2str(ort) '\circ)'],'t',[.08 .08 .84 .88]);
                 set(h,'FontSize',16)
                 switch saveFlag
                   case {1}
                     export_fig(pdfFilename,'-pdf','-nocrop','-append');
                 end

                %% Plots
                figure('name',['Bitmap generation: ' num2str(ort) char(0176)]) 

                % Plot columnar map and gaussian mask overlay
                subplot(2,5,1)
                [xContour,yContour]=find(~isnan(nanContour));
                bound=boundary(xContour,yContour);
                imgsc(bitmapCamspace);hold on; plot(yContour(bound),xContour(bound),'color','r','linewidth',.5);
                title({'Columns with contour overlay','(cam-space)'})
                % Plot targeted columns
                subplot(2,5,2); imgsc(bitmapCamspacePostMask);title({'Targeted columns','(cam-space)'})
                % plot resized
                subplot(2,5,3); imgsc(bitmapProjSpaceResized);title({'Imresize','(proj-space)'});  axis square; 
                % plot black-bars
                subplot(2,5,4); imgsc(bitmapProjSpaceBB);title({'Crop, black bars & flipud','(proj-space)'})
                % Plot post projector-camera alignment
                subplot(2,5,5); imgsc(bitmapProjSpaceAligned);title({'Camera-projector alignment','(proj-space)'})
                % Plot bitmap column area sizes
                subplot(2,5,6:10); 
                lenDC10=15;
                if ~isempty(columnAreas)
                    %bitmapStatsColumnDC=sqrt(columnAreas)/lenDC10;
                    %cdfplot(bitmapStatsColumnDC);grid off
                    xlabel('Estimated DC')
                    ylabel('Proportion %')
                    scatter(1:numel(columnAreas),sort(columnAreas,'descend'),50,'k','LineWidth',2)
                    legend({'BMP column area'})
                    title({['Column areas within contour' sprintf(' (n=%.0f)',numel(columnAreas))],...
                      sprintf('Total pixels (cam-space) = %.0f (%.1f%%)',...
                      bitmapData.pixelsON(bitmapNo),bitmapData.pixelsONDensity(bitmapNo))},...
                      'FontWeight','Normal')
                    xlabel('Columns')
                    ylabel('Area per column (pixels)')
                    xlim([0 numel(columnAreas)+1])
                    xticks([0:1:numel(columnAreas)+1])
                    legend({'Column'}); offwarning;
                    yticks([0:100:roundup(max(columnAreas),1000)])
                    ylim([0 roundup(max(columnAreas),1000)])
                end
                upFontSize(18,.0025)

           case {'proj2cam'} % Convert from camera space (512 x 512px) to projector space (1920 x 1080px)
                bitmapProjSpaceAligned=cloneBitmaps(currentBlockStruct,'load');

               %% Plot 1: Select central n-columns with a contour (from 2D gaussian)
                %[gaussianMask, bitmapCamspacePostMask] = isolateColumns(bitmapData, imagingData, blockID, bitmapNo, bitmapCamspace);
               
                %% Reverse camera-projector transformation
                bitmapProjSpaceUnaligned(:,:,bitmapNo)=applyCamProjAlignment(imagingData, blockID, bitmapProjSpaceAligned(:,:,bitmapNo+1), conversionType);

               %% Un-cropping (of y) and removing black bars (of x)
                bitmapProjSpaceUnBB(:,:,bitmapNo) = addBlackBars(bitmapData, bitmapProjSpaceUnaligned(:,:,bitmapNo), conversionType);

                %% Resize to cam space
                bitmapCamspace(:,:,bitmapNo) = resizeToSpace(imagingData, bitmapData, blockID, bitmapProjSpaceUnBB(:,:,bitmapNo), conversionType);

                %% Save to bitmapData
                % Calculate bitmap stats for columns within the contour (total pixels on, % pixel on, spatial duty cycle)
                contourROI=bitmapCamspace(:,:,bitmapNo);
               [labeledImage,columnAreas,boundaryBox]=getBlobAreas(contourROI);

                % save
                bitmapData.pixelsON(bitmapNo,blockID)=nansum(contourROI(:)>0);
                bitmapData.pixelsONDensity(bitmapNo,blockID)=nansum(contourROI(:)>0)*100/nansum(contourROI(:)==0); %prctPixON=((DC/(5/10))^2)/(2/100);
                bitmapData.columnarbitmapTFprojspace(:,:,bitmapNo,blockID)=bitmapProjSpaceAligned(:,:,bitmapNo+1); %first bitmap is black bmp
                bitmapData.columnarbitmapTFcamspace(:,:,bitmapNo,blockID)=bitmapCamspace(:,:,bitmapNo);
                bitmapData.nColumns(bitmapNo,blockID)=numel(columnAreas);
                bitmapData.medianColumnAreas(bitmapNo,blockID)=median(columnAreas);
                
                %% Power density calculations
                [bitmapData.adjustedSPD_uW(1,bitmapNo,blockID), bitmapData.ledpower_mW(1,bitmapNo,blockID)] = calculateSPD(behavioralData, imagingData, bitmapData, currentBlockStruct,...
                    bitmapNo,blockID,0);
                                
                %% Plots
                figure('name',['Bitmap generation: ' num2str(ort) char(0176)]) 

                % Plot columnar map and gaussian mask overlay
                subplot(2,5,1)
                %[xContour,yContour]=find(~isnan(nanContour));
                %bound=boundary(xContour,yContour);
                imgsc(bitmapCamspace(:,:,bitmapNo));hold on; %plot(yContour(bound),xContour(bound),'color','r','linewidth',.5);
                title({'Columns with contour overlay','(cam-space)'})
                % Plot targeted columns
                subplot(2,5,2); imgsc(bitmapCamspace(:,:,bitmapNo));title({'Targeted columns','(cam-space)'})
                % plot resized
                subplot(2,5,3); imgsc(bitmapCamspace(:,:,bitmapNo));title({'Imresize','(proj-space)'});  axis square; 
                % plot black-bars
                subplot(2,5,4); imgsc(bitmapProjSpaceUnBB(:,:,bitmapNo));title({'Crop, black bars & flipud','(proj-space)'})
                % Plot post projector-camera alignment
                subplot(2,5,5); imgsc(bitmapProjSpaceUnaligned(:,:,bitmapNo));title({'Camera-projector alignment','(proj-space)'})
                % Plot bitmap column area sizes
                subplot(2,5,6:10); 
                lenDC10=15;
                if ~isempty(columnAreas)
                    %bitmapStatsColumnDC=sqrt(columnAreas)/lenDC10;
                    %cdfplot(bitmapStatsColumnDC);grid off
                    xlabel('Estimated DC')
                    ylabel('Proportion %')
                    scatter(1:numel(columnAreas),sort(columnAreas,'descend'),50,'k','LineWidth',2)
                    legend({'BMP column area'})
                    title({['Column areas within contour' sprintf(' (n=%.0f)',numel(columnAreas))],...
                      sprintf('Total pixels (cam-space) = %.0f (%.1f%%)',...
                      bitmapData.pixelsON(bitmapNo),bitmapData.pixelsONDensity(bitmapNo))},...
                      'FontWeight','Normal')
                    xlabel('Columns')
                    ylabel('Area per column (pixels)')
                    xlim([0 numel(columnAreas)+1])
                    xticks([0:1:numel(columnAreas)+1])
                    legend({'Column'}); offwarning;
                    yticks([0:100:roundup(max(columnAreas),1000)])
                    ylim([0 roundup(max(columnAreas),1000)])
                end
                upFontSize(18,.0025)
    end
end
end
end













%% === SUBFUNCTIONS ===
%% Get masks based on gaussian-shaped responses to a visual gaussian. Contains n-contours that can be used to mask the original columnar map to isolate columns
function [imagingData, ort, bitmapCamspace, contourMaskLevel]=getColumnSelectionMask(bitmapData, imagingData, blockID, bitmapNo, plotFlag, saveFlag, pdfFilename)
    % Generate bitmap here

    % Extract necessary data from input structures
    ort = bitmapData.orts(1, bitmapNo, blockID);
    bitmapCamspace = bitmapData.columnarbitmapCoreg(:, :, bitmapNo, blockID); % Input image
    contourMaskLevel = bitmapData.gaussianContourLevel(1, bitmapNo, blockID); % Contour mask level

    % Fit a 2D Gaussian to create masks
    [imagingData] = fitGaussianResponseMask(bitmapData, imagingData, bitmapCamspace, blockID, bitmapNo, plotFlag);

    % Append to PDF if required
    appendToPDF(saveFlag, pdfFilename);
end

function [imagingData] = fitGaussianResponseMask(bitmapData, imagingData, bitmapCamspace, blockID, bitmapNo, plotFlag)
    if ~isnan(bitmapData.gaussianCond(blockID))
        gaussianResponse = abs(imagingData.gaussresp(:, :, bitmapData.gaussianCond(blockID), blockID));
        desiredResp = imresize(gaussianResponse, [imagingData.pixels(blockID), imagingData.pixels(blockID)], 'bilinear');
        [ROIMaskGaussian, centerCoords] = fit2Dgaussian(bitmapData.gaussianContourLevelMax(bitmapNo, blockID), desiredResp, bitmapCamspace, imagingData.nanmask(:, :, blockID), plotFlag);

        % Store fitted Gaussian mask and center coordinates
        imagingData.gaussfit(1, :, blockID) = padStruct(ROIMaskGaussian, 200);
        imagingData.centerCoords(1, :, blockID) = centerCoords;
    else
        imagingData.centerCoords(1, :, blockID) = [NaN, NaN]; % Default values if gaussianCond is NaN
    end
end

function appendToPDF(saveFlag, pdfFilename)
    switch saveFlag
        case 1
            export_fig(pdfFilename, '-pdf', '-nocrop', '-append');
    end
end

%% Isolate wanted columns using gaussian mask
function [gaussianMask, bitmapCamspacePostMask] = isolateColumns(bitmapData, imagingData, blockID, bitmapNo, bitmapCamspace)
    % Initialize outputs
    gaussianMask = [];
    bitmapCamspacePostMask = [];
    
    % Check if nColumnsWanted is not empty
    if ~isempty(bitmapData.nColumnsWanted)
        nColumnsWanted = bitmapData.nColumnsWanted(1, bitmapNo, blockID);
        contourMaskLevel = bitmapData.gaussianContourLevel(1, bitmapNo, blockID);
        disp(['Columns wanted: ' num2str(nColumnsWanted)]);
        
        % Single column, special case with contour level 200
        if nColumnsWanted == 1 && contourMaskLevel == 200
            disp('Loop 1');
            gaussianMask, bitmapCamspacePostMask = processSingleColumnSpecial(bitmapCamspace, imagingData.centerCoords(:, :, blockID));
        % Single or multiple columns, general case
        elseif (nColumnsWanted == 1 && contourMaskLevel < 200) || nColumnsWanted > 1
            disp('Loop 2');
            gaussianMask = imagingData.gaussfit(1, contourMaskLevel, blockID).area;
            bitmapCamspacePostMask = gaussianMask .* bitmapCamspace;
        
        % No columns wanted
        elseif nColumnsWanted == 0
            disp('Loop 3');
            gaussianMask = zeros(imagingData.pixels(1), imagingData.pixels(1));
            bitmapCamspacePostMask = gaussianMask .* bitmapCamspace;
        else
            disp('Loop 4');
            gaussianMask = ones(imagingData.pixels(1), imagingData.pixels(1));
            bitmapCamspacePostMask = gaussianMask .* bitmapCamspace;
        end
    end
end

function [gaussianMask, bitmapCamSpace] = processSingleColumnSpecial(inpict, centerCoords)
    % Process for a special single-column case with contour level 200
    L = bwlabel(inpict);
    S = regionprops(L, 'centroid');
    C = vertcat(S.Centroid); % Centers of all dots
    c0 = centerCoords; % Center of image
    DS = (C(:, 1) - c0(2)).^2 + (C(:, 2) - c0(1)).^2;
    [~, idx] = min(DS); % Minimize distance
    cp = C(idx, :); % This is the center point coordinates
    bitmapCamSpace = double(L == idx);
    gaussianMask = bitmapCamSpace;
end

%% Process image: Resize
function bitmapResized = resizeToSpace(imagingData, bitmapData, blockID, image,conversionType)
    % Resize to projector space
    % Calculate dimensions in projector space based on camera dimensions (in mm) and projector size (in px)
    switch conversionType
        case {'cam2proj'}
            targetPixelSize = imagingData.pixels(blockID)*imagingData.pixelsizemm(blockID) * 1 / bitmapData.pixelsizemm;
            bitmapResized = imresize(image, [targetPixelSize targetPixelSize], 'bilinear');
        case{'proj2cam'}
            targetPixelSize = imagingData.pixels(blockID);
            bitmapResized = imresize(image, [targetPixelSize targetPixelSize], 'bilinear');
    end
end

%% Process image: Add black bars (projector is wider than camera)
function image = addBlackBars(bitmapData, inputImage, conversionType)
switch conversionType
    case {'cam2proj'}
        % Calculate deltas for cropping (Y) and black bar sizes (X)
        imageProjSpace_deltaX = (bitmapData.projectorx - size(inputImage, 2)) / 2;
        imageProjSpace_deltaY = (bitmapData.projectory - size(inputImage, 1)) / 2;
    
        % Apply centered-crop the Y-dimension of the image if projector-Y < resizedCamImage
        if imageProjSpace_deltaY < 0
            imageProjSpace_preBB = inputImage(-imageProjSpace_deltaY + 1:end + imageProjSpace_deltaY, :);
        else
            imageProjSpace_preBB = inputImage;
        end
    
        % Add black bars to the X-dimension of the camera frame
        imageProjSpace_Y=size(imageProjSpace_preBB, 1);
        imageProjSpace_deltaX=round(abs(imageProjSpace_deltaX));
        blackbar_X = zeros(imageProjSpace_Y, imageProjSpace_deltaX);
        imageProjSpace_postBB = horzcat(blackbar_X, imageProjSpace_preBB, blackbar_X);
    
        % Ensure bitmap size is equal to projector (check for odd pixel)
        if mod(size(imageProjSpace_postBB, 2), 2) ~= 0
            imageProjSpace_postBB_rmpix = imageProjSpace_postBB(:, 1:end-1); % Ensure the image width matches the projector space width
        end
    
        % Optionally, adjust the bitmap for projector orientation
        imageProjSpace_postBB_flipud = flipud(imageProjSpace_postBB_rmpix);
        image=imageProjSpace_postBB_flipud;

    case {'proj2cam'} % recovers camera space image
        % Step 1: Flip the image back
        tmp=flipud(inputImage);
        % Step 2: Re-add the pixel removed from the end (tested, it is correct)
        tmp2=horzcat(tmp,zeros(bitmapData.projectory,1));

        % Step 3: Remove flanking bb=1080x199
        blackbars_x=1080;
        blackbars_y=199;
        tmp3=tmp2(:,blackbars_y+1:end-blackbars_y);
        % Step 4: Re-add black bars that were removed from first axis, and
        % remove pixel from the end (tested, it is correct)
        blackbars_y=zeros(round((size(tmp3,2)-size(tmp3,1)) / 2),size(tmp3,2));
        tmp4=vertcat(blackbars_y,tmp3,blackbars_y);
        tmp4cropped=tmp4(1:end-1,:);        
        image=tmp4cropped;
end
end

%% Apply cam-projector alignment
function bitmapProjSpaceAligned=applyCamProjAlignment(imagingData, blockID, bitmapProjSpaceBB, conversionType)
switch conversionType
    case {'cam2proj'}
            bitmapProjSpaceAligned = imwarp(bitmapProjSpaceBB,imagingData.transformmatrix(:,blockID),'OutputView',imref2d(size(bitmapProjSpaceBB)));
    case {'proj2cam'}
            bitmapProjSpaceAligned = imwarp(bitmapProjSpaceBB,invert(imagingData.transformmatrix(:,blockID)),'OutputView',imref2d(size(bitmapProjSpaceBB)));        
end
end

%% Calculate column areas
function [labeledImage,areasThresholded,boundaryBox]=getBlobAreas(binaryImage)
% define minimum area threshold (pesky dots)
lenDC01=15; %DC0.1 square
minAreaThreshold=lenDC01*2/3;
% label img blobs
labeledImage = bwlabel(binaryImage>0);%binaryImage);
% Measure the area
areas = struct2mat(regionprops(labeledImage, 'Area'));
boundaryBox = struct2mat(regionprops(labeledImage, 'BoundingBox'));
areasThresholded=areas(areas>minAreaThreshold); %remove pesky small dots
end