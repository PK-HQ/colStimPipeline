function [projspaceBitmapTform,bitmapStats,bitmapStats2]=convertForProjector(dataStructReference,dataStructSession,columnarBitmapCoregistered,Orts,...
    gridSize,gammaCorrFactor,sensitivity,conversionType,transformationMatrix,saveFlagBMP,saveFlag)

%% Load ROI mask
filenameStructSession=generateFilenames(dataStructSession);
filenameStructReference=generateFilenames(dataStructReference);
pdfFilename=filenameStructSession.neurometricPDF;
load(filenameStructReference.Orientation,'RespCond','Mask'); %contains mask, RespCondPCA

%% Paths
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
load(filenameStructReference.TS);%load([mainPath 'Chip/Chip20221104/run0/M28D20221104R0TS.mat'])

%% Set camera image-space and projector image-space dimensions
cam_size_px_CAMSPACE=TS.Header.Imaging.Resolution;
cam_px_size_mm=TS.Header.Imaging.SizePxl;
proj_px_x_PROJSPACE=1920;
proj_px_y_PROJSPACE=1080;
proj_size_px=.0054;
% Convert from cam-space pixel to mm space
cam_size_mm_x=cam_size_px_CAMSPACE*cam_px_size_mm; %512 * 0.0161
cam_size_mm_y=cam_size_mm_x;


% Set ROI mask zeros to nan
Mask=double(Mask);
ROIMaskNan=Mask; ROIMaskNan(ROIMaskNan==0)=NaN;

%% Generate bitmap here
for imgNo = 1:size(columnarBitmapCoregistered,3) % for each input image in camera space
    
    % Define image and contour map wanted
    ort=Orts(imgNo);
    camspaceBitmap=columnarBitmapCoregistered(:,:,imgNo); % grab input image wanted
    contourMaskLevel=dataStructSession.gaussianContourLevel(imgNo); % grab contour mask desired to extract n-columns
    
    %% Fit a 2D gaussian over gaussian response to create gaussian masks with an ROI mask overlay
    %if 2D gaussian already fitted, load. Else fit with data
    if isfile(filenameStructSession.gaussianFit) %load if gaussian fit exists
        % load transformaiton file if exist
        load(filenameStructSession.gaussianFit);
        if ~isequal(size(ROIMaskgaussian,2),contourMaskLevel) % if nContours is not same, fit
            % Load FFTed gaussian response
            load(dataStructSession.gaussianResponse)  
            desiredResp=imresize(abs(FFTCond(:,:,dataStructSession.gaussianCond)), [size(Mask,1) size(Mask,2)],'bilinear'); %FilterFermi2D(RespCond(:,:,1), desiredBandpass(1), desiredBandpass(2), TS.Header.Imaging.SizePxl);

            %fit and plot and save 2D gaussian with ROI mask overlay
            figure('name','Gaussian masks')
            [ROIMaskgaussian]=fit2Dgaussian(dataStructSession,desiredResp,camspaceBitmap,ROIMaskNan);
            save(filenameStructSession.gaussianFit,'ROIMaskgaussian')
            switch saveFlag
              case {1}
                export_fig(pdfFilename,'-pdf','-nocrop','-append');
            end
        end
    else % if no gaussian fit exists, fit
      if isfield(dataStructSession,'gaussianResponse') & isfile(dataStructSession.gaussianResponse)

        % Load FFTed gaussian response
        load(dataStructSession.gaussianResponse)  
        desiredResp=imresize(abs(FFTCond(:,:,dataStructSession.gaussianCond)), [size(Mask,1) size(Mask,2)],'bilinear'); %FilterFermi2D(RespCond(:,:,1), desiredBandpass(1), desiredBandpass(2), TS.Header.Imaging.SizePxl);

        %fit and plot and save 2D gaussian with ROI mask overlay
        figure('name','Gaussian masks')
        [ROIMaskgaussian]=fit2Dgaussian(dataStructSession,desiredResp,ROIMaskNan);
        save(filenameStructSession.gaussianFit,'ROIMaskgaussian')
        switch saveFlag
          case {1}
            export_fig(pdfFilename,'-pdf','-nocrop','-append');
        end
      else 
        % Dummy case
        ROIMaskgaussian(1).area=ones(size(columnarBitmapCoregistered(:,:,1),1), size(columnarBitmapCoregistered(:,:,1),2));
        dataStructSession.gaussianContourLevel=[1 1];
      end
    end
    
    if ort==0 || ort ==90 % Process only 0 or 90 deg orientations (for speed)
        switch conversionType
            case {'cam2proj'} % Convert from camera space (512 x 512px) to projector space (1920 x 1080px)
                
              %% Plot 1: Select central n-columns with a contour (from 2D gaussian)
                % Select gaussian mask, mask cam-space bitmap and restrict to n-columns
                gaussianMask = ROIMaskgaussian(contourMaskLevel).area;
                gaussianMaskedImg=gaussianMask .* camspaceBitmap;

                % Calculate bitmap stats for columns within the contour (total pixels on, % pixel on, spatial duty cycle)
                nanContour=double(gaussianMask);
                nanContour(nanContour==0)=NaN;
                contourROI=nanContour.*camspaceBitmap;
                bitmapStats.contourPix(imgNo)=nansum(contourROI(:)>0);
                bitmapStats.contourPixDensity(imgNo)=nansum(contourROI(:)>0)*100/nansum(contourROI(:)==0); %prctPixON=((DC/(5/10))^2)/(2/100);
                bitmapStats.contourDC(imgNo)=sqrt(bitmapStats.contourPixDensity(imgNo)*2/100)*5/10; %DC=sqrt(prctPixON*2/100)*5/10;
                
                % Same but for the entire ROI mask (d'>6)
                bitmapStats.maskPixDensity(imgNo)=nansum(contourROI(:)>0)*100/nansum(ROIMaskNan(:));
                bitmapStats.maskDC(imgNo)=sqrt(bitmapStats.maskPixDensity(imgNo)*2/100)*5/10;

                figure('name',['Projector alignment: ' num2str(ort) char(0176) ' bitmapStats']) 
                subplot(2,6,1)
                [xContour,yContour]=find(~isnan(nanContour));
                bound=boundary(xContour,yContour);
                imgsc(camspaceBitmap);hold on; plot(yContour(bound),xContour(bound),'color','r','linewidth',.5);
                title({'Columns with contour overlay','(cam-space)'})

                subplot(2,6,2)
                imgsc(gaussianMaskedImg);title({'Targeted columns','(cam-space)'})


               %% Plot 2: Resize to projector space
                %project to projector space
                cam_size_px_x_PROJSPACE=cam_size_mm_x*1/proj_size_px;
                cam_size_px_y_PROJSPACE=cam_size_mm_y*1/proj_size_px;
                img_px_PROJSPACE=imresize(gaussianMaskedImg, [cam_size_px_x_PROJSPACE cam_size_px_y_PROJSPACE],'bilinear');
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
                projspaceBitmapTform=projspaceBitmapAdjusted;
                subplot(2,6,5)
                imgsc(projspaceBitmapAdjusted);title({'Crop, black bars & flipud','(proj-space)'})

                
                %% Apply camera-projector transformation
                if ~isempty(transformationMatrix)
                    projspaceBitmapTform = imwarp(projspaceBitmapAdjusted,transformationMatrix,'OutputView',imref2d(size(img_px_BB_PROJSPACE)));
                end

                subplot(2,6,6)
                imgsc(projspaceBitmapTform);title({'Camera-projector alignment','(proj-space)'})

                %{
                %resize
                x_cent = round(size(columnarBitmapCoregisteredAdjusted,1)/2);
                y_cent = round(size(columnarBitmapCoregisteredAdjusted,2)/2);
                sizeX = size(img_px_BB_PROJSPACE,1);
                sizeY = size(img_px_BB_PROJSPACE,2);
                %centroide = [x_cent y_cent]; % this isn't used
                %I2 = imcrop(I,rect) crops the image I. rect is a four-element position vector of the
                %form [xmin ymin width height] that specifies the size and position of the crop rectangle. 
                %imcrop returns the cropped image, I2.
                xmin = x_cent-sizeX/2;
                ymin = y_cent-sizeY/2;
                columnarBitmapCoregisteredTform = imcrop(columnarBitmapCoregisteredAdjusted,[xmin+1 ymin+1 sizeX-1 sizeY-1]);
                %}
                if size(projspaceBitmapTform,1)==1080 & size(projspaceBitmapTform,2)==1920
                    %
                else
                    fprintf('Warning: Check bitmapStats size, not 1080x1920\n')
                end


                %% Check size relative to known DCs

                subplot(2,6,7:12)
                binaryImage=projspaceBitmapTform;%bitmapStatsRef(:,:,1);

                if size(binaryImage,1)==512
                    proj2camFactor=.0054/.0161;
                else
                    proj2camFactor=1;
                end


                % convert these to DC equivalent, then plot hist with DC size x Count
                % Label the image so we can get the average perpendicular width.
                lenDC10=15;
                [labeledImage,bitmapStats2,boundaryBox]=getBlobAreas(binaryImage,lenDC10*2/3);
                bitmapStats.nBlobs(imgNo)=numel(bitmapStats2);
                bitmapStats.medianBlobAreas(imgNo)=median(bitmapStats2);

                %% Plot column area and estimated DC
                bitmapStatsColumnDC=sqrt(bitmapStats2)/lenDC10;
                cdfplot(bitmapStatsColumnDC);grid off
                xlabel('Estimated DC')
                ylabel('Proportion %')
                scatter(1:numel(bitmapStats2),sort(bitmapStats2,'descend'),'k','filled')
                DC10area=lenDC10*lenDC10;
                DC20area=(lenDC10*2)^2;
                DC30area=(lenDC10*3)^2;
                DC15area=(lenDC10*1.5)^2;
                modulationSuccess=2500;
                line(xlim,[DC10area DC10area]*proj2camFactor,'color',.7*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[DC20area DC20area]*proj2camFactor,'color',.5*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[DC30area DC30area]*proj2camFactor,'color',.3*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                %line(xlim,[modulationSuccess modulationSuccess]*proj2camFactor,'color',[46,139,87]/255,'lineWidth',1.5,'LineStyle','--')
                legend({'BMP column area', 'DC0.1 area', 'DC0.2 area', 'DC0.3 area', 'Good modulations'})
                title({['Columns by area ' sprintf('(n=%.0f)',numel(bitmapStats2))],...
                  sprintf('[Contour area] Pixels=%.0f (%.1f%%), DC_{spatial}=%.2f',...
                  bitmapStats.contourPix(imgNo),bitmapStats.contourPixDensity(imgNo),bitmapStats.contourDC(imgNo)),...
                  sprintf('[ROI mask area] Pixels=%.0f (%.1f%%), DC_{spatial}=%.2f',...
                  bitmapStats.contourPix(imgNo),bitmapStats.maskPixDensity(imgNo),bitmapStats.maskDC(imgNo))})
                xlabel('Columns')
                ylabel('Pixels-on per column')
                xlim([1 numel(bitmapStats2)+1])
                xticks([1:1:numel(bitmapStats2)+1])
                h = drawrectangle('Position',[0,DC10area,numel(bitmapStats2),DC30area-DC10area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                legend({'Column','DC0.1','DC0.2','DC0.3','DC0.success','','','',}); offwarning;
                upFontSize(14,.0025)
                yticks([0:1000:roundup(max(bitmapStats2),1000)])
                ylim([0 roundup(max(bitmapStats2),1000)])

                %% Save columnar optostim BMP
                if ispc & saveFlagBMP==1

                    bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                    bmpFilename = sprintf('O%05gHE%04gG%03gS%05gC%02g.bmp',...
                        ort*100,gridSize,gammaCorrFactor,sensitivity*10000,dataStructSession.gaussianContourLevel(imgNo));

                    if ~exist(bmpPath, 'dir')
                         mkdir(bmpPath)
                    end

                    % Save the bitmapStats
                    disp(bmpFilename)
                    imwrite(projspaceBitmapTform,[bmpPath bmpFilename]);
                end

                %% Save the figures
                 [~,h]=suplabel(['Converting bitmap for projector dimensions (' num2str(ort) '\circ)'],'t',[.08 .08 .84 .88]);
                 set(h,'FontSize',16)
                 switch saveFlag
                   case {1}
                     export_fig(pdfFilename,'-pdf','-nocrop','-append');
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
                    fprintf('Warning: Check bitmapStats size, not 1080x1920\n')
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
                bitmapStats2=struct2mat(measurements);
                bitmapStats2=bitmapStats2(bitmapStats2>50); %remove pesky small dots
                scatter(1:numel(bitmapStats2),sort(bitmapStats2,'descend'),'k','filled')
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
                title(sprintf('Sorted bitmapStats column areas (Pix %=%.2g, DC=%.2g, DCt=)',bitmapStats.pixelDensity(imgNo),bitmapStats.estDutycycle(imgNo)))
                xlabel('Sorted blobs')
                ylabel('Pixel-ON area')
                xlim([0 numel(bitmapStats2)])
                h = drawrectangle('Position',[0,DC15area,numel(bitmapStats2),2500-DC15area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                upFontSize(14,.01)

                %% Save columnar optostim BMP
                bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                bmpFilename = sprintf('O%05gHE%04gG%03gT%05g.bmp',...
                    ort*100,gridSize,gammaCorrFactor,sensitivity*100);
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
                camspaceBitmap=imtranslate(camspaceBitmap,round([-projX -projY]))

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

%{



% 8.22mm x  8.22mm
cam_mm_x=CamPxX*CamPxSize;
cam_mm_y=CamPxY*CamPxSize;
% 10.368mm x 5.832mm
proj_mm_x=ProjPxX*ProjPxSize;
proj_mm_y=ProjPxY*ProjPxSize;

%Camera in projector pixel space (1522 x 1522)
cam_in_proj_space_px_x=cam_mm_x/ProjPxSize; %projector x pixels that camera covers
cam_in_proj_space_px_y=cam_mm_y/ProjPxSize; %projector y pixels that camera covers (includes the vertical black bars)

%Projector in camera pixel space (645 x 363)
proj_in_cam_space_px_x=proj_mm_x/CamPxSize; %camera x pixels that projector covers
proj_in_cam_space_px_y=proj_mm_y/CamPxSize; %camera y pixels that projector covers

switch conversionType
    case {'cam2proj'}
       %% Load files
        filenameStruct=generateFilenames(dataStructSession);
        load(filenameStruct.TS);


        camImg=imread(filenameStruct.calImg);
        projBMP=imread(filenameStruct.calBmp);
        if size(camImg,3)>1
            camImg=rgb2gray(camImg);
        end
        if size(projBMP,3)>1
            projBMP=rgb2gray(projBMP);
        end
        %% Convert bmp from camera px to projector px

        %% Figure: Overlay of projector and camera windows
        figure; axis equal
        xlim([0 ProjPxX])
        ylim([0 ProjPxX])
        drawrectangle('Position',[(ProjPxX-cam_in_proj_space_px_x)/2,...
            abs(cam_in_proj_space_px_y-ProjPxX)/2,cam_in_proj_space_px_x,cam_in_proj_space_px_y],...
            'Label','Camera Area','Color','k');
        drawrectangle('Position',[0,abs(ProjPxY-ProjPxX)/2,ProjPxX,ProjPxY],...
            'Label','Projector Area','Color','r');
        xlabel('X-pixels')
        ylabel('Y-pixels')
        set(gcf,'color','w');
        title('Overlay of projector and camera windows')

        %% Figure: Resize columnar map from Cam-PX to Proj-PX
        colBMPProjPx = imresize(inputImg,[cam_in_proj_space_px_x cam_in_proj_space_px_y]); %CamBMP in mm space
        figure; axis equal
        xlim([0 ProjPxX])
        ylim([0 ProjPxX])
        imshowpair(colBMPProjPx,projBMP)
        xlabel('X-pixels')
        ylabel('Y-pixels')
        set(gcf,'color','w');

        %% Find xy translate and rotation of projector relative to camera frame
        % Get centers of images
        szCalImg=size(colBMPProjPx);
        midCamFrameX=(szCalImg(1)+1)/2;
        midCamFrameY=(szCalImg(2)+1)/2;
        midProjFrameX=795; % predefined with selection %256.5;
        midProjFrameY=739; %238;

        % Apply xy translate, rotate (leftward shift = -x, cw rotate = -angle)
        proj_px_offset_x=midCamFrameX-midProjFrameX;%midCamFrameX-midProjFrameX;
        proj_px_offset_y=midCamFrameY-midProjFrameY; 
        projRot=-2.5;%-2.5;

        colBMPProjPxR=imrotate(colBMPProjPx,projRot);
        %[midProjFrameX,midProjFrameY,~] = impixelzoom(colBMPProjPxR);


        colBMPProjPxTR=imtranslate(colBMPProjPxR,[proj_px_offset_x, proj_px_offset_y]);
        subplot(1,3,1)
        imagesc(colBMPProjPx);axis square;
        line([0 szCalImg(1)], [szCalImg(1)/2 szCalImg(1)/2],'color','r','lineWidth',1)
        line([szCalImg(1)/2 szCalImg(1)/2], [0 szCalImg(1)],'color','r','lineWidth',1)
        title('Coregistered projector to camera frame')
        subplot(1,3,2)
        imagesc(colBMPProjPxR);axis square;
        line([0 szCalImg(1)], [szCalImg(1)/2 szCalImg(1)/2],'color','r','lineWidth',1)
        line([szCalImg(1)/2 szCalImg(1)/2], [0 szCalImg(1)],'color','r','lineWidth',1)
        title('Translate')
        subplot(1,3,3)
        imagesc(colBMPProjPxTR);axis square;
        line([0 szCalImg(1)], [szCalImg(1)/2 szCalImg(1)/2],'color','r','lineWidth',1)
        line([szCalImg(1)/2 szCalImg(1)/2], [0 szCalImg(1)],'color','r','lineWidth',1)
        title('Translate + Rotate')
        suplabel(['bitmapStats for ' num2str(ort)],'t')
        %Crop Y add X
        cam_remove_dim1_length=(size(colBMPProjPxTR,1)-size(projBMP,1))/2;
        cam_add_dim2_length=(size(colBMPProjPxTR,2)-size(projBMP,2))/2;

        if size(colBMPProjPxTR,3)>1
            colBMPProjPxTR=rgb2gray(colBMPProjPxTR);
        end
        projbitmapStatsTRBB=colBMPProjPxTR;
        projbitmapStatsTRBB(1:cam_remove_dim1_length,:)=[];
        projbitmapStatsTRBB(end-cam_remove_dim1_length+1:end,:)=[];
        projbitmapStatsTRBB=horzcat(zeros(abs(cam_add_dim2_length),ProjPxY)', projbitmapStatsTRBB);
        projbitmapStatsTRBB=horzcat(projbitmapStatsTRBB,zeros(abs(cam_add_dim2_length),ProjPxY)');

        %% flip ud because projector is upsidedown
        projbitmapStatsTRBB=flipud(projbitmapStatsTRBB);
        projbitmapStatsTRBB_unflipped=(projbitmapStatsTRBB);
        imshowpair(projbitmapStatsTRBB,projBMP);
    case {'proj2cam'}
        %reverse translate and rotate
        
        %clip horizontal BB
        
        %add vertical bb
       
end
 %}

end