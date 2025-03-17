% Load, extract O0 and O90 for contrast zero, then bandpass filter for columns
monkeyName='Chip';
monkeyID='28';
date='20240118';
ortmapRun='0';
optoRun='10';
optofixRun='5';

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

% Load, extract O0 and O90 for contrast zero, then bandpass filter for columns
load(['Y:\' monkeyName '\' monkeyName  date '\run' optofixRun '\M' monkeyID 'D' date 'R' optofixRun 'StabIntgS004E010.mat'])
load(['Y:\' monkeyName '\' monkeyName  date '\run' optofixRun '\M' monkeyID 'D' date 'R' optofixRun 'TS.mat'])
optoConds=[];
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
O0=find(optoConds==0); %select 0
O90=find(optoConds==90); %select 90
% get images
O0img=columnarFilter(TS,DataCond(:,:,O0)) .* nanMask; %filter
O90img=columnarFilter(TS,DataCond(:,:,O90)) .* nanMask; %filter
superImg=cat(3, O0Vximg, O90Vximg, O0img, O90img);
[mse0, ssimValue0, ncc0] = imageSimilarity(O0Vximg, O0img);
[mse90, ssimValue90, ncc90] = imageSimilarity(O90Vximg, O90img);

figure('Name','State similarity')
nRows=2;nCols=3;
gap=.15;marginV=.01;marginH=.075;
[hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.13 marginV+.13], [marginH marginH]);
axes(hAx(1));
imgsc(O0Vximg); axis square; title('Task-state 0\circ'); colorbar; cax(.9, superImg); upFontSize(20,.015); addPix2MM(O90Vximg,1,2,3) %colormap(hAx(1), fireice); 
axes(hAx(2)); 
imgsc(O0img); axis square; title('Fix-state 0\circ');  colorbar; cax(.9, superImg); upFontSize(20,.015); addPix2MM(O90Vximg,2,2,3)  %colormap(hAx(2), fireice); 
axes(hAx(3)); 
imgsc(double(imfuse(O0Vximg,O0img, 'method', 'blend'))/255 .* nanMask); axis square; title({'Task-Fix 0\circ', sprintf('NCC: %.2f',ncc0)}); hColorbar=colorbar; ylabel(hColorbar, 'Difference (Normalized)'); upFontSize(20,.015); addPix2MM(O90Vximg,3,2,3) 
upFontSize(24,.015)
axes(hAx(4));
imgsc(O90Vximg); axis square; title('Task-state 90\circ'); colorbar; cax(.9, superImg); upFontSize(20,.015); addPix2MM(O90Vximg,4,2,3)  %colormap(hAx(4),fireice); 
axes(hAx(5)); 
imgsc(O90img); axis square; title('Fix-state 90\circ'); colorbar; cax(.9, superImg); upFontSize(20,.015); addPix2MM(O90Vximg,5,2,3)  %colormap(hAx(5),fireice);
axes(hAx(6)); 
imgsc(double(imfuse(O90Vximg,O90img, 'method', 'blend'))/255 .* nanMask); axis square; title({'Task-Fix 90\circ', sprintf('NCC: %.2f',ncc90)}); hColorbar=colorbar; ylabel(hColorbar, 'Difference (Normalized)'); upFontSize(20,.015); addPix2MM(O90Vximg,6,2,3); upFontSize(20,.015);

% Set all font properties to ensure consistency
set(findall(gcf, '-property', 'FontUnits'), 'FontUnits', 'points');
set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'Tex');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Arial');
export_fig(['Y:\' monkeyName '\Meta\stateSimilarity\M28D' date 'R' optoRun 'R' optofixRun], '-pdf', '-append', '-nocrop', '-painters');
export_fig(['Y:\' monkeyName '\Meta\stateSimilarity\M28D' date 'R' optoRun 'R' optofixRun], '-svg', '-append', '-nocrop', '-painters');
export_fig(['Y:\' monkeyName '\Meta\stateSimilarity\M28D' date 'R' optoRun 'R' optofixRun], '-png', '-append', '-nocrop', '-painters');