entryNo=4;
dataStruct(entryNo).monkeyNo='28';
dataStruct(entryNo).monkey='Chip';
dataStruct(entryNo).date='20220920';
dataStruct(entryNo).run='0';
dataStruct(entryNo).site='30';
dataStruct(entryNo).modality='GCaMP';

filenameStruct=generateFilenames(dataStruct(4));
load(filenameStruct.TS);

%% Def size params
cam_size_px_CAMSPACE=TS.Header.Imaging.Resolution;
cam_px_size_mm=TS.Header.Imaging.SizePxl;
proj_px_x=1920;
proj_px_y=1080;
proj_size_px=.0054;

% Convert camera image from pixel to mm space
cam_size_mm_x=cam_size_px_CAMSPACE*cam_px_size_mm; %512 * 0.0161
cam_size_mm_y=cam_size_mm_x;

%project to projector space
cam_in_proj_space_px_x=cam_size_mm_x*1/proj_size_px;
cam_in_proj_space_px_y=cam_size_mm_y*1/proj_size_px;

%calculate projector misalignment
cam_px_mid_x=(cam_in_proj_space_px_x+1)/2;
cam_px_mid_y=(cam_in_proj_space_px_y+1)/2;
proj_px_offset_x=795; % predefined with selection %256.5;
proj_px_offset_y=739; %238;

%% Figure: Overlay of projector and camera windows
figure;

% True xy translate, rotate (leftward shift = -x, cw rotate = -angle)
projRotTrue=-2.5;
proj_px_delta_x_True=cam_px_mid_x-proj_px_offset_x;
proj_px_delta_y_True=(cam_px_mid_y-proj_px_offset_y); 

subplot(1, 3, 1)
axis equal
xlim([0 proj_px_x])
ylim([0 proj_px_x])
drawrectangle('Position',[(proj_px_x-cam_in_proj_space_px_x)/2,...
    abs(cam_in_proj_space_px_y-proj_px_x)/2,cam_in_proj_space_px_x,cam_in_proj_space_px_y],...
    'Label','Camera Area','Color','k');
drawrectangle('Position',[0+proj_px_delta_x_True,(abs(proj_px_y-proj_px_x)/2)+proj_px_delta_y_True,proj_px_x,proj_px_y],...
    'Label','Projector Area','Rotation',projRotTrue,'Color','r');
xlabel('X-pixels')
ylabel('Y-pixels')
set(gcf,'color','w');
title('True projector-camera alignment')

%% Corrected
% True xy translate, rotate (leftward shift = -x, cw rotate = -angle)
projRotCorrected=projRotTrue-(projRotTrue);
proj_px_delta_x_Corrected=proj_px_delta_x_True-(proj_px_delta_x_True);
proj_px_delta_y_Corrected=proj_px_delta_y_True-(proj_px_delta_y_True); 

subplot(1, 3, 2)
axis equal
xlim([0 proj_px_x])
ylim([0 proj_px_x])
drawrectangle('Position',[(proj_px_x-cam_in_proj_space_px_x)/2,...
    abs(cam_in_proj_space_px_y-proj_px_x)/2,cam_in_proj_space_px_x,cam_in_proj_space_px_y],...
    'Label','Camera Area','Color','k');
drawrectangle('Position',[0+proj_px_delta_x_Corrected,(abs(proj_px_y-proj_px_x)/2)+proj_px_delta_y_Corrected,proj_px_x,proj_px_y],...
    'Label','Projector Area','Rotation',projRotCorrected,'Color','r');
xlabel('X-pixels')
ylabel('Y-pixels')
set(gcf,'color','w');
title(sprintf('Corrected projector-camera alignment (Rot=%0.3g,X=%0.3g,Y=%0.3g)',(projRotTrue),(proj_px_delta_x_True),(proj_px_delta_y_True)))

%% Test image
subplot(1, 3, 3)
axis equal
xlim([0 proj_px_x])
ylim([0 proj_px_x])
xlabel('X-pixels')
ylabel('Y-pixels')
set(gcf,'color','w');
title(sprintf('Superimposed pre-post corrected camera image')
