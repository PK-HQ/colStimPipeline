% Load, extract O0 and O90 for contrast zero, then bandpass filter for columns
monkeyName='Chip';
monkeyID='28';
date='20230829';
ortmapRun='0';
optoRun='3';

load(['Y:\' monkeyName '\' monkeyName date '\run' ortmapRun '\M' monkeyID 'D' date 'R' ortmapRun 'OrientationP2.mat'])
load(['Y:\' monkeyName '\' monkeyName date '\run' optoRun '\M' monkeyID 'D' date 'R' optoRun 'StabIntgS004E010.mat'])
load(['Y:\' monkeyName '\' monkeyName date '\run' optoRun '\M' monkeyID 'D' date 'R' optoRun 'TS.mat'])

% Make d' mask
nanMask = double(Mask); nanMask(~Mask) = NaN; % Replace 0s with NaN

gaborCon=TS.Header.Conditions.StimCon;
gaborOrt=TS.Header.Conditions.GaborOrt;
for i = 1:numel(TS.Header.Conditions.ProjImg)
    if contains(TS.Header.Conditions.ProjImg{i}, 'O000')
        optoConds(i) = 0;
    elseif contains(TS.Header.Conditions.ProjImg{i}, 'O090')
        optoConds(i) = 90;
    else
        optoConds(i) = NaN;
    end
end

%get conditions
O0Vx=find(gaborCon==0 & optoConds==0); %select 0
O90Vx=find(gaborCon==0 & optoConds==90); %select 90
% get images
O0Vximg=columnarFilter(TS,mean(DataCond(:,:,O0Vx),3)) .* nanMask; %filter
O90Vximg=columnarFilter(TS,mean(DataCond(:,:,O90Vx),3)) .* nanMask; %filter

superImg=cat(3, O0Vximg, O90Vximg);

figure('Name','State similarity')
nRows=1;nCols=2;
gap=.15;marginV=.01;marginH=.075;
[hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.13 marginV+.13], [marginH marginH]);
axes(hAx(1));
imgsc(O0Vximg); axis square; title('Task-state 0\circ'); colorbar; cax(.4, superImg); upFontSize(20,.015); addPix2MM(O90Vximg,1,1,2) %colormap(hAx(1), fireice); 
upFontSize(24,.015)
axes(hAx(2));
imgsc(O90Vximg); axis square; title('Task-state 90\circ'); colorbar; cax(.4, superImg); upFontSize(20,.015); addPix2MM(O90Vximg,2,1,2)  %colormap(hAx(4),fireice); 
upFontSize(24,.015)

% Set all font properties to ensure consistency
set(findall(gcf, '-property', 'FontUnits'), 'FontUnits', 'points');
set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'Tex');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
export_fig(['Y:\' monkeyName '\Meta\neurometric\M28D' date 'R' optoRun '-novisual'], '-pdf', '-append', '-nocrop', '-painters');
export_fig(['Y:\' monkeyName '\Meta\neurometric\M28D' date 'R' optoRun '-novisual'], '-svg', '-append', '-nocrop', '-painters');
export_fig(['Y:\' monkeyName '\Meta\stateSimilarity\M28D' date 'R' optoRun '-novisual'], '-png', '-append', '-nocrop', '-painters');