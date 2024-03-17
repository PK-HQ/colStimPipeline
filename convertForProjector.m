function [bitmapData]=convertForProjector(behavioralData, imagingData, bitmapData,...
    currentBlockStruct, conversionType, blockID, ...
    pdfFilename, plotFlag,saveFlagBMP,saveFlag)
disp('Generating bitmap...')
%% Set camera image-space and projector image-space dimensions
cam_size_px_CAMSPACE=imagingData.pixels(blockID);%behavioralData.Header.Imaging.Resolution;
cam_px_size_mm=imagingData.pixelsizemm(blockID); %behavioralData.Header.Imaging.SizePxl;
proj_px_x_PROJSPACE=bitmapData.projectorx;
proj_px_y_PROJSPACE=bitmapData.projectory;
proj_size_px=bitmapData.pixelsizemm;
% Convert from cam-space pixel to mm space
cam_size_mm_x=cam_size_px_CAMSPACE*cam_px_size_mm; %512 * 0.0161
cam_size_mm_y=cam_size_mm_x;

%% Generate bitmap here
for bitmapNo = 1:size(bitmapData.columnarbitmapCoreg,3) % for each input image in camera space
    
    % Define image and contour map wanted
    ort=bitmapData.orts(1,bitmapNo,blockID);
    camspaceBitmap=bitmapData.columnarbitmapCoreg(:,:,bitmapNo,blockID); % grab input image wanted
    contourMaskLevel=bitmapData.gaussianContourLevel(1,bitmapNo,blockID); % grab contour mask desired to extract n-columns
    
    %% Fit a 2D gaussian over gaussian response to create gaussian masks with an ROI mask overlay
        % fit
        if ~isnan(bitmapData.gaussianCond(blockID))
            gaussianResponse=abs(imagingData.gaussresp(:,:,bitmapData.gaussianCond(blockID),blockID));
            desiredResp=imresize(gaussianResponse,...
                [imagingData.pixels(blockID) imagingData.pixels(blockID)],...
                'bilinear');
            [ROIMaskgaussian,centerCoords]=fit2Dgaussian(bitmapData.gaussianContourLevelMax(bitmapNo,blockID),desiredResp,camspaceBitmap,imagingData.nanmask(:,:,blockID),plotFlag);
            %save(currentBlockStruct.gaussianFit,'ROIMaskgaussian','centerCoords')
            imagingData.gaussfit(1,:,blockID)=padStruct(ROIMaskgaussian, 200);
            imagingData.centerCoords(1,:,blockID)=centerCoords;
        else
            %imagingData.gaussfit(1,:,blockID)=repmat(NaN, 1, 200);
            %imagingData.centerCoords(1,:,blockID)=[NaN NaN];
        end
        
        switch saveFlag
          case {1}
            export_fig(pdfFilename,'-pdf','-nocrop','-append');
        end
        %isempty(imagingData.gaussfit(1,contourMaskLevel,blockID)) || isnan(imagingData.gaussfit(1,contourMaskLevel,blockID).threshPercentile) || isnan(imagingData.centerCoords(:,blockID))
        %end

    if ort==0 || ort ==90 % Process only 0 or 90 deg orientations (for speed)
        switch conversionType
            case {'cam2proj'} % Convert from camera space (512 x 512px) to projector space (1920 x 1080px)
                
                
               %% Plot 1: Select central n-columns with a contour (from 2D gaussian)
                if ~isempty(bitmapData.nColumnsWanted)
                    disp(['Columns wanted:' num2str(bitmapData.nColumnsWanted(1,bitmapNo,blockID))])
                    if bitmapData.nColumnsWanted(1,bitmapNo,blockID)==1 && bitmapData.gaussianContourLevel(1,bitmapNo,blockID)==200 % Single column 
                        disp('Loop 1')
                        % new test
                        inpict = camspaceBitmap;
                        L = bwlabel(inpict);
                        S = regionprops(L,'centroid');
                        C = vertcat(S.Centroid); % centers of all dots
                        c0 = (imagingData.centerCoords); % center of image
                        % squared euclidean distance from all dots to center
                        DS = (C(:,1)-c0(2)).^2 + (C(:,2)-c0(1)).^2;
                        [~,idx] = min(DS); % minimize
                        % this is the center point coordinates and image
                        cp = C(idx,:);
                        bitmapCamSpace=double(L==idx);
                        gaussianMask = bitmapCamSpace;
                    elseif bitmapData.nColumnsWanted(1,bitmapNo,blockID)==1 && bitmapData.gaussianContourLevel(1,bitmapNo,blockID)<200  % single-column
                        disp('Loop 2')
                        % Select gaussian mask, mask cam-space bitmap and restrict to n-columns
                        gaussianMask = imagingData.gaussfit(1,contourMaskLevel,blockID).area;
                        bitmapCamSpace=gaussianMask .* camspaceBitmap;
                    elseif bitmapData.nColumnsWanted(1,bitmapNo,blockID)>1  % multi-column
                        disp('Loop 3')
                        % Select gaussian mask, mask cam-space bitmap and restrict to n-columns
                        gaussianMask = imagingData.gaussfit(1,contourMaskLevel,blockID).area;
                        bitmapCamSpace=gaussianMask .* camspaceBitmap;
                    elseif bitmapData.nColumnsWanted(1,bitmapNo,blockID)==0
                        disp('Loop 4')
                        gaussianMask = zeros(imagingData.pixels(1),imagingData.pixels(1));
                        bitmapCamSpace=gaussianMask .* camspaceBitmap;
                    elseif isnan(imagingData.centerCoords(1,1,blockID)) || isnan(bitmapData.nColumnsWanted(1,bitmapNo,blockID))
                        disp('Loop 5')
                        gaussianMask = ones(imagingData.pixels(1),imagingData.pixels(1));
                        bitmapCamSpace=gaussianMask .* camspaceBitmap;
                    end
                end

                % Calculate bitmap stats for columns within the contour (total pixels on, % pixel on, spatial duty cycle)
                nanContour=double(gaussianMask);
                nanContour(nanContour==0)=NaN;
                contourROI=nanContour.*camspaceBitmap;
                bitmapData.pixelsON(bitmapNo,blockID)=nansum(contourROI(:)>0);
                bitmapData.pixelsONDensity(bitmapNo,blockID)=nansum(contourROI(:)>0)*100/nansum(contourROI(:)==0); %prctPixON=((DC/(5/10))^2)/(2/100);
                %bitmapData.pixelsONDutyCycle(1,bitmapNo,blockID)=sqrt(bitmapData.pixelsONDensity(bitmapNo)*2/100)*5/10; %DC=sqrt(prctPixON*2/100)*5/10;
                
                % Same but for the entire ROI mask (d'>6)
                %bitmapData.maskPixDensity(1,bitmapNo,blockID)=nansum(contourROI(:)>0)*100/nansum(imagingData.nanmask(:));
                %bitmapData.maskDC(1,bitmapNo,blockID)=sqrt(bitmapData.maskPixDensity(bitmapNo)*2/100)*5/10;

                figure('name',['Bitmap generation: ' num2str(ort) char(0176)]) 
                subplot(2,6,1)
                [xContour,yContour]=find(~isnan(nanContour));
                bound=boundary(xContour,yContour);
                imgsc(camspaceBitmap);hold on; plot(yContour(bound),xContour(bound),'color','r','linewidth',.5);
                if exist('c0','var')
                    plot(c0(2),c0(1),'+','DisplayName','Gaussian center');
                    plot(cp(1),cp(2),'.','DisplayName','Nearest column');
                end
                title({'Columns with contour overlay','(cam-space)'})

                % Plot targeted columns
                subplot(2,6,2)
                imgsc(bitmapCamSpace);title({'Targeted columns','(cam-space)'})
                
               %% Plot 2: Resize to projector space
                %project to projector space
                cam_size_px_x_PROJSPACE=cam_size_mm_x*1/proj_size_px;
                cam_size_px_y_PROJSPACE=cam_size_mm_y*1/proj_size_px;
                img_px_PROJSPACE=imresize(bitmapCamSpace, [cam_size_px_x_PROJSPACE cam_size_px_y_PROJSPACE],'bilinear');
                subplot(2,6,3)
                imgsc(img_px_PROJSPACE);title({'Imresize','(proj-space)'})

              %% Correction for rotate and translate
                %calculate projector misalignment
                cam_px_mid_x=(cam_size_px_x_PROJSPACE+1)/2;
                cam_px_mid_y=(cam_size_px_y_PROJSPACE+1)/2;
                proj_px_offset_x=795;
                proj_px_offset_y=739;

                % Apply xy translate, rotate (leftward shift = -x, cw rotate = -angle)
                calibrationRotate=0;
                calibrationX=0;
                calibrationY=0;

                projRot=0-calibrationRotate;
                proj_px_delta_x=(cam_px_mid_x-proj_px_offset_x)+calibrationX;
                proj_px_delta_y=(cam_px_mid_y-proj_px_offset_y)-calibrationY; 

                %correct for projector rotation(CW=-angle) > translation (left = -x) > flipud
                img_px_R_PROJSPACE=imrotate(img_px_PROJSPACE,projRot,'crop','bilinear'); %Linear interpolation, supposed to be better than nearest neighbour
                img_px_RT_PROJSPACE=imtranslate(img_px_R_PROJSPACE,[proj_px_delta_x, proj_px_delta_y],...
                    'bilinear',... %Linear interpolation, supposed to be better than nearest neighbour
                    'FillValues',0,... %fills gaps with 0, crops image to same as original size
                    'OutputView','same');
                subplot(2,6,4)
                imgsc(img_px_RT_PROJSPACE);title({'Rotate & Translate','(inter-space)'})

                %convert to grayscale if rgb
                if size(img_px_RT_PROJSPACE,3)>1
                    img_px_RT_PROJSPACE=rgb2gray(img_px_RT_PROJSPACE);
                end

               %% Plot 3: Cropping (of y) and adding black bars (of x)
                %crop the y of camera frame
                cam_x_delta=(proj_px_x_PROJSPACE-size(img_px_RT_PROJSPACE,1))/2;
                cam_y_delta=-(proj_px_y_PROJSPACE-size(img_px_RT_PROJSPACE,2))/2;
                img_px_BB_PROJSPACE=img_px_RT_PROJSPACE(cam_y_delta+1:end-cam_y_delta,:);

                % add black bars to x of camera frame
                blackbar=zeros(size(img_px_BB_PROJSPACE,1),round(cam_x_delta));
                img_px_BB_PROJSPACE=horzcat(blackbar,img_px_BB_PROJSPACE);
                img_px_BB_PROJSPACE=horzcat(img_px_BB_PROJSPACE,blackbar);
                img_px_BB_PROJSPACE=img_px_BB_PROJSPACE(:,1:end-1); %crop to 1920, ACCOUNT FOR WITH A TRANSLATION OR REMOVE ODD NUMBER AFTER RESIZE!
                projspaceBitmapAdjusted=flipud(img_px_BB_PROJSPACE);
                bitmapProjSpaceTransformed=projspaceBitmapAdjusted;
                subplot(2,6,5)
                imgsc(projspaceBitmapAdjusted);title({'Crop, black bars & flipud','(proj-space)'})

                %% Apply camera-projector transformation
                if ~isempty(imagingData.transformmatrix(:,blockID))
                    bitmapProjSpaceTransformed = imwarp(projspaceBitmapAdjusted,imagingData.transformmatrix(:,blockID),'OutputView',imref2d(size(img_px_BB_PROJSPACE)));
                end

                subplot(2,6,6)
                imgsc(bitmapProjSpaceTransformed);title({'Camera-projector alignment','(proj-space)'})

                %% Check size relative to known DCs

                subplot(2,6,7:12)
                binaryImage=bitmapProjSpaceTransformed;%bitmapStatsRef(:,:,1);

                if size(binaryImage,1)==imagingData.pixels
                    proj2camFactor=bitmapData.pixelsizemm/imagingData.pixelsizemm;
                else
                    proj2camFactor=1;
                end


                % convert these to DC equivalent, then plot hist with DC size x Count
                % Label the image so we can get the average perpendicular width.
                lenDC10=15;
                [labeledImage,columnAreas,boundaryBox]=getBlobAreas(binaryImage,lenDC10*2/3);
                bitmapData.nColumns(bitmapNo,blockID)=numel(columnAreas);
                bitmapData.medianColumnAreas(bitmapNo,blockID)=median(columnAreas);

                %% Plot column area and estimated DC
                if ~isempty(columnAreas)
                    bitmapStatsColumnDC=sqrt(columnAreas)/lenDC10;
                    cdfplot(bitmapStatsColumnDC);grid off
                    xlabel('Estimated DC')
                    ylabel('Proportion %')
                    scatter(1:numel(columnAreas),sort(columnAreas,'descend'),'k','filled')
                    DC10area=lenDC10*lenDC10;
                    DC20area=(lenDC10*2)^2;
                    DC30area=(lenDC10*3)^2;
                    DC15area=(lenDC10*1.5)^2;
                    modulationSuccess=2500;
                    %line(xlim,[DC10area DC10area]*proj2camFactor,'color',.7*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                    %line(xlim,[DC20area DC20area]*proj2camFactor,'color',.5*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                    %line(xlim,[DC30area DC30area]*proj2camFactor,'color',.3*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                    %line(xlim,[modulationSuccess modulationSuccess]*proj2camFactor,'color',[46,139,87]/255,'lineWidth',1.5,'LineStyle','--')
                    legend({'BMP column area', 'DC0.1 area', 'DC0.2 area', 'DC0.3 area', 'Good modulations'})
                    title({['Columns by area ' sprintf('(n=%.0f)',numel(columnAreas))],...
                      sprintf('ROI area, Pixels-ON=%.0f (%.1f%%)',...
                      bitmapData.pixelsON(bitmapNo),bitmapData.pixelsONDensity(bitmapNo))})
                    xlabel('Columns')
                    ylabel('Area per column (pixels)')
                    xlim([0 numel(columnAreas)+1])
                    xticks([0:1:numel(columnAreas)+1])
                    %h = drawrectangle('Position',[0,DC10area,numel(columnAreas),DC30area-DC10area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                    legend({'Column','DC0.1','DC0.2','DC0.3','DC0.success','','','',}); offwarning;
                    upFontSize(14,.0025)
                    yticks([0:1000:roundup(max(columnAreas),1000)])
                    ylim([0 roundup(max(columnAreas),1000)])
                end

                %% Save columnar optostim BMP
                bitmapData.columnarbitmapTFprojspace(:,:,bitmapNo,blockID)=bitmapProjSpaceTransformed;
                bitmapData.columnarbitmapTFcamspace(:,:,bitmapNo,blockID)=bitmapCamSpace;
                if ispc & saveFlagBMP==1
                    bmpPath=sprintf('T:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                    bmpFilename = sprintf('O%05gHE%04gG%03gS%05gC%02g.bmp',...
                        ort*100,bitmapData.gridSize,bitmapData.gammaCorrFactor,bitmapData.sensitivity*10000,bitmapData.gaussianContourLevel(bitmapNo));

                    if ~exist(bmpPath, 'dir')
                         mkdir(bmpPath)
                    end

                    % Save the bitmapData
                    disp(bmpFilename)
                    imwrite(bitmapProjSpaceTransformed,[bmpPath bmpFilename]);
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
        end
    end
end
                 
                 
                 
                 
                 
            case {'calibrate'}
                %imread(filenameStruct.bmpOriginal);
                %imread(filenameStruct.responseOriginal);

                % Convert camera image from pixel to mm space
                cam_size_mm_x=cam_size_px_CAMSPACE*cam_px_size_mm; %512 * 0.0161
                cam_size_mm_y=cam_size_mm_x;

                %project to projector space
                cam_size_px_x_PROJSPACE=cam_size_mm_x*1/proj_size_px;
                cam_size_px_y_PROJSPACE=cam_size_mm_y*1/proj_size_px;
                img_px_PROJSPACE=imresize(camspaceBitmap, [cam_size_px_x_PROJSPACE cam_size_px_y_PROJSPACE],'bilinear');

                %calculate projector misalignment
                cam_px_mid_x=(cam_size_px_x_PROJSPACE+1)/2;
                cam_px_mid_y=(cam_size_px_y_PROJSPACE+1)/2;
                proj_px_offset_x=795; % predefined with selection %256.5;
                proj_px_offset_y=739; %238;

                % Apply xy translate, rotate (leftward shift = -x, cw rotate = -angle)

                calibrationRotate=-4.1-.461;
                calibrationX=27-3.9;
                calibrationY=-35-17.7;

                projRot=-2.5-calibrationRotate; %minus the rotation value from calibration
                %orig
                proj_px_delta_x=(cam_px_mid_x-proj_px_offset_x)+calibrationX;
                proj_px_delta_y=(cam_px_mid_y-proj_px_offset_y)-calibrationY; 

                %% Correction for rotate and translate
                %correct for projector rotation(CW=-angle) > translation (left = -x) > flipud
                img_px_R_PROJSPACE=imrotate(img_px_PROJSPACE,projRot,'crop','bilinear'); %Linear interpolation, supposed to be better than nearest neighbour
                img_px_RT_PROJSPACE=imtranslate(img_px_R_PROJSPACE,[proj_px_delta_x, proj_px_delta_y],...
                    'bilinear',... %Linear interpolation, supposed to be better than nearest neighbour
                    'FillValues',0,'OutputView','same'); %fills gaps with 0, crops image to same as original size

                %convert to grayscale if rgb
                if size(img_px_RT_PROJSPACE,3)>1
                    img_px_RT_PROJSPACE=rgb2gray(img_px_RT_PROJSPACE);
                end

                %% Cropping (y) and adding black bars (x)
                %crop the y of camera frame
                cam_x_delta=(proj_px_x_PROJSPACE-size(img_px_RT_PROJSPACE,1))/2;
                cam_y_delta=-(proj_px_y_PROJSPACE-size(img_px_RT_PROJSPACE,2))/2;
                img_px_BB_PROJSPACE=img_px_RT_PROJSPACE(cam_y_delta+1:end-cam_y_delta,:);

                % add black bars to x of camera frame
                blackbar=zeros(size(img_px_BB_PROJSPACE,1),round(cam_x_delta));
                img_px_BB_PROJSPACE=horzcat(blackbar,img_px_BB_PROJSPACE);
                img_px_BB_PROJSPACE=horzcat(img_px_BB_PROJSPACE,blackbar);
                img_px_BB_PROJSPACE=img_px_BB_PROJSPACE(:,1:end-1); %crop to 1920, ACCOUNT FOR WITH A TRANSLATION OR REMOVE ODD NUMBER AFTER RESIZE!
                projspaceBitmapAdjusted=flipud(img_px_BB_PROJSPACE);

                if size(projspaceBitmapAdjusted,1)==1080 & size(projspaceBitmapAdjusted,2)==1920
                    %
                else
                    fprintf('Warning: Check bitmapData size, not 1080x1920\n')
                end


                %% Check size relative to known DCs
                figure
                binaryImage=projspaceBitmapAdjusted;%bitmapStatsRef(:,:,1);

                if size(binaryImage,1)==512
                    proj2camFactor=.0054/.0161;
                else
                    proj2camFactor=1;
                end
                % Label the image so we can get the average perpendicular width.
                labeledImage = bwlabel(binaryImage);
                % Measure the area
                measurements = regionprops(labeledImage, 'Area');
                columnAreas=struct2mat(measurements);
                columnAreas=columnAreas(columnAreas>50); %remove pesky small dots
                scatter(1:numel(columnAreas),sort(columnAreas,'descend'),'k','filled')
                DC10area=15*15;
                DC20area=(15*2)^2;
                DC30area=(15*3)^2;
                DC15area=(15*1.5)^2;
                modulationSuccess=2500;
                line(xlim,[DC10area DC10area]*proj2camFactor,'color',.7*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[DC20area DC20area]*proj2camFactor,'color',.5*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[DC30area DC30area]*proj2camFactor,'color',.3*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[modulationSuccess modulationSuccess]*proj2camFactor,'color',[46,139,87]/255,'lineWidth',1.5,'LineStyle','--')
                legend({'BMP column area', 'DC0.1 area', 'DC0.2 area', 'DC0.3 area', 'Good modulations'})
                title(sprintf('Sorted bitmapData column areas (Pix %=%.2g, DC=%.2g, DCt=)',bitmapData.pixelDensity(bitmapNo),bitmapData.estDutycycle(bitmapNo)))
                xlabel('Sorted blobs')
                ylabel('Pixel-ON area')
                xlim([0 numel(columnAreas)])
                h = drawrectangle('Position',[0,DC15area,numel(columnAreas),2500-DC15area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                upFontSize(14,.01)

                %% Save columnar optostim BMP
                bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                bmpFilename = sprintf('O%05gHE%04gG%03gT%05g.bmp',...
                    ort*100,bitmapData.gridSize,bitmapData.gammaCorrFactor,bitmapData.sensitivity*100);
                %exportgraphics(projbitmapStatsTRBB,bmpFilename);%,'Resolution',300)

                if ~exist(bmpPath, 'dir')
                     mkdir(bmpPath)
                end

                switch saveFlag
                    case {1}
                       %imwrite(camspaceBitmapAdjusted,[bmpPath bmpFilename]);
                    case {0}
                end
            case {'proj2cam'}
                %flipud
                camspaceBitmap=flipud(camspaceBitmap);

                %-X +Y

                %translate
                camspaceBitmap=imtranslate(camspaceBitmap,round([-projX -projY]));

                %rotate
                camspaceBitmap=imrotate(camspaceBitmap, -projRot);

                %*projpxsz
                [x, y]=size(camspaceBitmap);
                camX=x*proj_px_sz/cam_px_sz;
                camY=y*proj_px_sz/cam_px_sz;

                camspaceBitmap=imresize(camspaceBitmap, [camX camY]);

        end
    end
end
if plotFlag==0
    close
end
end