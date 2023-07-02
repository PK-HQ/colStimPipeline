function [columnarBitmapsTform,nBlobs,medianBlobAreas]=convertForProjector(dataStruct,columnarBitmaps,Orts,...
    gridSize,gammaCorrFactor,sensitivity,conversionType,transformationMatrix,saveOrNot)


filenameStruct=generateFilenames(dataStruct);
pdfFilename=filenameStruct.neurometricPDF;

if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
load([mainPath 'Chip/Chip20221104/run0/M28D20221104R0TS.mat'])
%load(filenameStruct.TS);


% Dims
cam_size_px_CAMSPACE=TS.Header.Imaging.Resolution;
cam_px_size_mm=TS.Header.Imaging.SizePxl;
proj_px_x_PROJSPACE=1920;
proj_px_y_PROJSPACE=1080;
proj_size_px=.0054;
for imgNo = 1:size(columnarBitmaps,3)
    inputImg=columnarBitmaps(:,:,imgNo);
    ort=Orts(imgNo);
    if ort==0 || ort ==90
        switch conversionType
            case {'cam2proj'}
                % Convert camera image from pixel to mm space
                cam_size_mm_x=cam_size_px_CAMSPACE*cam_px_size_mm; %512 * 0.0161
                cam_size_mm_y=cam_size_mm_x;


                figure('name',['Projector alignment:' num2str(ort) ' bitmap']) 
                subplot(2,5,1)
                imgsc(inputImg);title('imresize')


                %project to projector space
                cam_size_px_x_PROJSPACE=cam_size_mm_x*1/proj_size_px;
                cam_size_px_y_PROJSPACE=cam_size_mm_y*1/proj_size_px;
                img_px_PROJSPACE=imresize(inputImg, [cam_size_px_x_PROJSPACE cam_size_px_y_PROJSPACE],'bilinear');


                subplot(2,5,2)
                imgsc(img_px_PROJSPACE);title('imresize')

                %calculate projector misalignment
                cam_px_mid_x=(cam_size_px_x_PROJSPACE+1)/2;
                cam_px_mid_y=(cam_size_px_y_PROJSPACE+1)/2;
                proj_px_offset_x=795; % predefined with selection %256.5;
                proj_px_offset_y=739; %238;

                % Apply xy translate, rotate (leftward shift = -x, cw rotate = -angle)
                calibrationRotate=0;%0-1-.1%-3.86*1.065; %0.7264 points per deg
                calibrationX=0;%55-17+5%78.2*.55; %0.0690
                calibrationY=0;%30-67+5%-96.8*.55; 

                projRot=0-calibrationRotate;%-2.5-calibrationRotate; %minus the rotation value from calibration
                %orig
                proj_px_delta_x=(cam_px_mid_x-proj_px_offset_x)+calibrationX;
                proj_px_delta_y=(cam_px_mid_y-proj_px_offset_y)-calibrationY; 

                %% Correction for rotate and translate
                %correct for projector rotation(CW=-angle) > translation (left = -x) > flipud
                img_px_R_PROJSPACE=imrotate(img_px_PROJSPACE,projRot,'crop','bilinear'); %Linear interpolation, supposed to be better than nearest neighbour
                img_px_RT_PROJSPACE=imtranslate(img_px_R_PROJSPACE,[proj_px_delta_x, proj_px_delta_y],...
                    'bilinear',... %Linear interpolation, supposed to be better than nearest neighbour
                    'FillValues',0,'OutputView','same'); %fills gaps with 0, crops image to same as original size


                subplot(2,5,3)
                imgsc(img_px_RT_PROJSPACE);title('R,T')


                %load('D:/OIData/Calibration20221025/transformParams.mat')
                %img_px_RT_PROJSPACE = imwarp(img_px_RT_PROJSPACE,transformParams,'OutputView',imref2d(size(img_px_RT_PROJSPACE)));

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
                columnarBitmapsAdjusted=flipud(img_px_BB_PROJSPACE);
                columnarBitmapsTform=columnarBitmapsAdjusted;

                subplot(2,5,4)
                imgsc(columnarBitmapsAdjusted);title('Crop+BB+flip')

                %apply transform for session
                if ~isempty(transformationMatrix)
                    columnarBitmapsTform = imwarp(columnarBitmapsAdjusted,transformationMatrix,'OutputView',imref2d(size(img_px_BB_PROJSPACE)));
                end

                subplot(2,5,5)
                imgsc(columnarBitmapsTform);title('Tform')

                %{
                %resize
                x_cent = round(size(columnarBitmapsAdjusted,1)/2);
                y_cent = round(size(columnarBitmapsAdjusted,2)/2);
                sizeX = size(img_px_BB_PROJSPACE,1);
                sizeY = size(img_px_BB_PROJSPACE,2);
                %centroide = [x_cent y_cent]; % this isn't used
                %I2 = imcrop(I,rect) crops the image I. rect is a four-element position vector of the
                %form [xmin ymin width height] that specifies the size and position of the crop rectangle. 
                %imcrop returns the cropped image, I2.
                xmin = x_cent-sizeX/2;
                ymin = y_cent-sizeY/2;
                columnarBitmapsTform = imcrop(columnarBitmapsAdjusted,[xmin+1 ymin+1 sizeX-1 sizeY-1]);
                %}
                if size(columnarBitmapsTform,1)==1080 & size(columnarBitmapsTform,2)==1920
                    %
                else
                    fprintf('Warning: Check bitmap size, not 1080x1920\n')
                end


                %% Check size relative to known DCs

                subplot(2,5,6:10)
                binaryImage=columnarBitmapsTform;%bitmapRef(:,:,1);

                if size(binaryImage,1)==512
                    proj2camFactor=.0054/.0161;
                else
                    proj2camFactor=1;
                end


                % convert these to DC equivalent, then plot hist with DC size x Count
                % Label the image so we can get the average perpendicular width.
                lenDC10=15;
                [labeledImage,bitmapColumnAreas]=getBlobAreas(binaryImage,lenDC10*2/3);
                nBlobs(imgNo)=numel(bitmapColumnAreas);
                medianBlobAreas(imgNo)=median(bitmapColumnAreas);
                
                %{
                labeledImage = bwlabel(binaryImage);
                % Measure the area
                measurements = regionprops(labeledImage, 'Area');
                bitmapColumnAreas=struct2mat(measurements);
                bitmapColumnAreas=bitmapColumnAreas(bitmapColumnAreas>10); %remove pesky small dots
                %}
                
                bitmapColumnDC=sqrt(bitmapColumnAreas)/lenDC10;
                cdfplot(bitmapColumnDC);grid off
                xlabel('Estimated DC')
                ylabel('Proportion %')

                scatter(1:numel(bitmapColumnAreas),sort(bitmapColumnAreas,'descend'),'k','filled')
                DC10area=lenDC10*lenDC10;
                DC20area=(lenDC10*2)^2;
                DC30area=(lenDC10*3)^2;
                DC15area=(lenDC10*1.5)^2;
                modulationSuccess=2500;
                line(xlim,[DC10area DC10area]*proj2camFactor,'color',.7*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[DC20area DC20area]*proj2camFactor,'color',.5*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[DC30area DC30area]*proj2camFactor,'color',.3*[1 1 1],'lineWidth',1.5,'LineStyle','--')
                line(xlim,[modulationSuccess modulationSuccess]*proj2camFactor,'color',[46,139,87]/255,'lineWidth',1.5,'LineStyle','--')
                legend({'BMP column area', 'DC0.1 area', 'DC0.2 area', 'DC0.3 area', 'Good modulations'})
                title(sprintf('Sorted bitmap column areas (n-blobs=%.0f)',numel(bitmapColumnAreas)))
                xlabel('Sorted blobs')
                ylabel('Pixel-ON area')
                xlim([0 numel(bitmapColumnAreas)])
                %h = drawrectangle('Position',[0,DC15area,numel(bitmapColumnAreas),2500-DC15area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                h = drawrectangle('Position',[0,DC15area,numel(bitmapColumnAreas),DC30area-DC15area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                legend({'','','','','','','','',})
                upFontSize(14,.01)
                ylim([0 7000])

                %% Save columnar optostim BMP
                if ispc
                    
                    bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                    bmpFilename = sprintf('O%05gHE%04gG%03gS%05g.bmp',...
                        ort*100,gridSize,gammaCorrFactor,sensitivity*10000);
                    
                    if ~exist(bmpPath, 'dir')
                         mkdir(bmpPath)
                    end
                        
                    % Save the bitmap
                    imwrite(columnarBitmapsTform,[bmpPath bmpFilename]);
                end

                %% Save the figures
                [~,h]=suplabel(['Converting bitmap for projector dimensions (' num2str(ort) '\circ)'],'t',[.08 .08 .84 .87]);
                set(h,'FontSize',16)
                export_fig(pdfFilename,'-pdf','-nocrop','-append');

            case {'calibrate'}
                %imread(filenameStruct.bmpOriginal);
                %imread(filenameStruct.responseOriginal);

                % Convert camera image from pixel to mm space
                cam_size_mm_x=cam_size_px_CAMSPACE*cam_px_size_mm; %512 * 0.0161
                cam_size_mm_y=cam_size_mm_x;

                %project to projector space
                cam_size_px_x_PROJSPACE=cam_size_mm_x*1/proj_size_px;
                cam_size_px_y_PROJSPACE=cam_size_mm_y*1/proj_size_px;
                img_px_PROJSPACE=imresize(inputImg, [cam_size_px_x_PROJSPACE cam_size_px_y_PROJSPACE],'bilinear');

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
                columnarBitmapsAdjusted=flipud(img_px_BB_PROJSPACE);

                if size(columnarBitmapsAdjusted,1)==1080 & size(columnarBitmapsAdjusted,2)==1920
                    %
                else
                    fprintf('Warning: Check bitmap size, not 1080x1920\n')
                end


                %% Check size relative to known DCs
                figure
                binaryImage=columnarBitmapsAdjusted;%bitmapRef(:,:,1);

                if size(binaryImage,1)==512
                    proj2camFactor=.0054/.0161;
                else
                    proj2camFactor=1;
                end
                % Label the image so we can get the average perpendicular width.
                labeledImage = bwlabel(binaryImage);
                % Measure the area
                measurements = regionprops(labeledImage, 'Area');
                bitmapColumnAreas=struct2mat(measurements);
                bitmapColumnAreas=bitmapColumnAreas(bitmapColumnAreas>10); %remove pesky small dots
                scatter(1:numel(bitmapColumnAreas),sort(bitmapColumnAreas,'descend'),'k','filled')
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
                title('Sorted bitmap column areas')
                xlabel('Sorted blobs')
                ylabel('Pixel-ON area')
                xlim([0 numel(bitmapColumnAreas)])
                h = drawrectangle('Position',[0,DC15area,numel(bitmapColumnAreas),2500-DC15area],'Color',[46,139,87]/255,'Label','Goldilocks zone');
                upFontSize(14,.01)

                %% Save columnar optostim BMP
                bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
                bmpFilename = sprintf('O%05gHE%04gG%03gT%05g.bmp',...
                    ort*100,gridSize,gammaCorrFactor,sensitivity*100);
                %exportgraphics(projBitmapTRBB,bmpFilename);%,'Resolution',300)

                if ~exist(bmpPath, 'dir')
                     mkdir(bmpPath)
                end

                switch saveOrNot
                    case {1}
                       %imwrite(columnarBitmapsAdjusted,[bmpPath bmpFilename]);
                    case {0}
                end
            case {'proj2cam'}
                %flipud
                inputImg=flipud(inputImg);

                %-X +Y

                %translate
                inputImg=imtranslate(inputImg,round([-projX -projY]))

                %rotate
                inputImg=imrotate(inputImg, -projRot);

                %*projpxsz
                [x, y]=size(inputImg);
                camX=x*proj_px_sz/cam_px_sz;
                camY=y*proj_px_sz/cam_px_sz;

                inputImg=imresize(inputImg, [camX camY]);

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
        filenameStruct=generateFilenames(dataStruct);
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
        suplabel(['Bitmap for ' num2str(ort)],'t')
        %Crop Y add X
        cam_remove_dim1_length=(size(colBMPProjPxTR,1)-size(projBMP,1))/2;
        cam_add_dim2_length=(size(colBMPProjPxTR,2)-size(projBMP,2))/2;

        if size(colBMPProjPxTR,3)>1
            colBMPProjPxTR=rgb2gray(colBMPProjPxTR);
        end
        projBitmapTRBB=colBMPProjPxTR;
        projBitmapTRBB(1:cam_remove_dim1_length,:)=[];
        projBitmapTRBB(end-cam_remove_dim1_length+1:end,:)=[];
        projBitmapTRBB=horzcat(zeros(abs(cam_add_dim2_length),ProjPxY)', projBitmapTRBB);
        projBitmapTRBB=horzcat(projBitmapTRBB,zeros(abs(cam_add_dim2_length),ProjPxY)');

        %% flip ud because projector is upsidedown
        projBitmapTRBB=flipud(projBitmapTRBB);
        projBitmapTRBB_unflipped=(projBitmapTRBB);
        imshowpair(projBitmapTRBB,projBMP);
    case {'proj2cam'}
        %reverse translate and rotate
        
        %clip horizontal BB
        
        %add vertical bb
       
end
 %}

end