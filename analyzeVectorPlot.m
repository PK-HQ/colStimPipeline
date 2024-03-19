%% Reverse mapAmpOrt, split into 12 ort maps by preferred ort
% Image and cond params
Height=size(MapOrt,1);
Width=size(MapOrt,2);
OrtRange=(Ort(2)-Ort(1))/2;
mapOrtSplit=nan(Height,Width,numel(Ort));
ROINaNMask=double(Mask); ROINaNMask(ROINaNMask==0)=NaN;
MapOrt=MapOrt.*ROINaNMask;
% Scale by amp
tAmpOrt = AmpOrt/max(AmpOrt(:))*3;
tAmpOrt(tAmpOrt>1) = 1; % why clip like that? 
tAmpOrt(MapOrt==-10) = 1; % why assign like that?

% Visualize
nRows=3;
nCols=4;
figure('name','12 orts')
[hAx,~]=tight_subplot(nRows,nCols,[.025 .025]);
for ort = 1:numel(Ort)
    ortColumns = double(MapOrt>=Ort(ort)-OrtRange & MapOrt<=Ort(ort)+OrtRange);
    ortColumns(ortColumns==0) = NaN;
    mapOrtSplit(:,:,ort)=ortColumns.*tAmpOrt;
    axes(hAx(ort));
    imgsc(mapOrtSplit(:,:,ort)); colorbar;
    addPix2MM(1,512,1,512,ort,nRows,nCols)
end
[~,h]=suplabel('12 orientations','t',[.08 .08 .84 .76]); set(h,'FontWeight','normal');

%% Load ort map, but separated into each color
 [angle, amplitude, angleAverage, amplitudeAverage]=calculateProjectionVector(condResp, mapOrtSplit, Ort);

% Plot circle map (one circle map per condition type: con opto, incon opto,
% visual)
figure
[figure_handle,count,speeds,directions,Table,Others] = WindRose(angles,amplitude);

plot(vector(:,1),vector(:,2))
plot(vectorAverage(1),vectorAverage(2))


% Overlay vector sum amplitude
