load('Y:\Chip\Chip20230803\run0\M28D20230803R0OrientationP2.mat','RespCond','RespCondPCA','Mask')
load('Y:\Chip\Chip20230803\run0\M28D20230803R0TS.mat','TS')

figure
[hAx,~]=tight_subplot(3,4,[.1 .15],[.05 .15],[.04 .04]);
% Raw responses
axes(hAx(1))
imgsc(RespCond(:,:,1),'Raw response');
axes(hAx(2))
imgsc(FilterFermi2D(RespCond(:,:,1), 0, .8, TS.Header.Imaging.SizePxl),'Low-pass: 0.0-0.8cpmm');
axes(hAx(3))
imgsc(FilterFermi2D(RespCond(:,:,1), 0.8, 3.0, TS.Header.Imaging.SizePxl),'Band-pass');

% Desired response for fitting
Mask=double(Mask);
ROIMaskNan=Mask; ROIMaskNan(ROIMaskNan==0)=NaN;
desiredBandpass=[0 .2];
desiredResp=FilterFermi2D(RespCond(:,:,1), desiredBandpass(1), desiredBandpass(2), TS.Header.Imaging.SizePxl);
desiredRespMask=desiredResp.*Mask;

axes(hAx(4))
imgsc(desiredRespMask,sprintf('Bandpass = %.1f to %.1f cpmm',desiredBandpass(1),desiredBandpass(2)));

% Fit 2d gaussian
fit2Dgaussian(desiredResp,ROIMaskNan)













%% Fit 2D gaussian
opts.iso=false;
opts.tilted=true;

% fit for bandpass 0-0.3cpmm
[xCoords,yCoords] = meshgrid(1:512);%meshgrid(-10:10,-20:20);
zAmplitudes = desiredRespMask; %exp(-(xi-3).^2-(yi+4).^2) + randn(size(xi));
results = autoGaussianSurf(xCoords,yCoords,zAmplitudes,opts);
% plot 2d gaussian
axes(hAx(7))
imgsc(flipud(results.G),'Estimated Gaussian: 0-0.3cpmm') % flip because of matlab's imagesc oddity

% fit for bandpass 0-0.5cpmm
[xCoords,yCoords] = meshgrid(1:512);%meshgrid(-10:10,-20:20);
zAmplitudes = FilterFermi2D(RespCond(:,:,1), 0, .4, TS.Header.Imaging.SizePxl); %exp(-(xi-3).^2-(yi+4).^2) + randn(size(xi));
results = autoGaussianSurf(xCoords,yCoords,zAmplitudes,opts);
% plot 2d gaussian
axes(hAx(8))
imgsc(flipud(results.G),'Estimated Gaussian: 0-0.4cpmm') % flip because of matlab's imagesc oddity

% fit for bandpass 0-0.79cpmm
[xCoords,yCoords] = meshgrid(1:512);%meshgrid(-10:10,-20:20);
zAmplitudes = FilterFermi2D(RespCond(:,:,1), 0, .7, TS.Header.Imaging.SizePxl); %exp(-(xi-3).^2-(yi+4).^2) + randn(size(xi));
results = autoGaussianSurf(xCoords,yCoords,zAmplitudes,opts);
% plot 2d gaussian
axes(hAx(9))
imgsc(flipud(results.G),'Estimated Gaussian: 0-0.7cpmm') % flip because of matlab's imagesc oddity