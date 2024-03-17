% Modified from Gero Nootz (2012), MATLAB FEX
%
% Coeficients A convention:
%	A = [Amplitude, x0, x-Width, y0, y-Width, Angle(in Radians)]
%
% X-data convention:
%	X is of size(n,n,2) where 
%	X(:,:,1) : x-coordinates,  
%	X(:,:,2) : y-coordinates.
%



function [ROIMaskgaussian,centerCoords]=fit2Dgaussian(varargin)
figure('name','Column selection')
nContourLevels=varargin{1};
desiredResp=varargin{2};
if length(varargin)==5
    nContourLevels=varargin{1};
    desiredResp=varargin{2};
    bitmap=varargin{3};
    mask=varargin{4};
    plotFlag=varargin{5};
end
%filenameStructSession=generateFilenames(dataStruct);
%desiredResp=transformImage(desiredResp,filenameStructSession);
% In this numerical test we use two-dimensional fitting functions:
% 1. 2D Gaussian function ( A requires 5 coefs ).
gaussFunc = @(A,X) A(1)*exp( -((X(:,:,1)-A(2)).^2/(2*A(3)^2) + (X(:,:,2)-A(4)).^2/(2*A(5)^2)) );
% 2. 2D Rotated Gaussian function ( A requires 6 coefs ).
rotgaussFunc = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );
%% ---Parameters---
n = size(desiredResp,1); m = size(desiredResp,2);         % n x m pixels area/data matrix
initialParams = [1,200,100,250,100,0.7854];   % Inital (guess) parameters
InterpMethod='nearest'; % 'nearest','linear','spline','cubic'
FitOrientation='fit';	% 'fit': fit for orientation, 'dont' fit for orientation
%% ---Build numerical Grids---
% Numerical Grid
[x,y]=meshgrid(1:n,1:m); X=zeros(m,n,2); X(:,:,1)=x; X(:,:,2)=y;
% High Resolution Grid
h=1; [xh,yh]=meshgrid(1:1/h:n,1:1/h:m); Xh=zeros(h*m,h*n,2); Xh(:,:,1)=xh; Xh(:,:,2)=yh;

%% ---Fit---
% Define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
lb = [0,1,0,1,0,0];
ub = [realmax('double'),n,(n)^2,n,(n)^2,pi/4];
% Fit sample data
%opts = optimset('Display','off');
opts = optimoptions('lsqcurvefit','Display','off');

switch FitOrientation
    case 'dont' 
        [A,resnorm,res,flag,output] = lsqcurvefit(gaussFunc,initialParams(1:5),X,desiredResp,lb(1:5),ub(1:5),opts);
    case 'fit'
        [A,resnorm,res,flag,output] = lsqcurvefit(rotgaussFunc,initialParams,X,desiredResp,lb,ub,opts);
    otherwise
        error('invalid entry');
end
%disp(output); % display summary of LSQ algorithm
%% ---Plot Data---
% Plot 3D Data and Fitted curve
%hf1=figure('name','3D gaussian'); set(hf1,'Position',[1000 600 800 500]); 
switch FitOrientation
    case 'dont' 
        C=del2(gaussFunc(A,X)); mesh(xh,yh,gaussFunc(A,X),C); hold on
        fullGaussianFit=gaussFunc(A,X);
    case 'fit'  
        C=del2(rotgaussFunc(A,X)); mesh(xh,yh,rotgaussFunc(A,X),C); hold on
        fullGaussianFit=rotgaussFunc(A,X);
end

%% Plot actual data
subplot(3,8,[1 2 9 10]); 
imagesc(x(1,:),y(:,1),desiredResp); axis square; 
colormap('gray');colorbar; 
title('4Hz response to gaussian');
upFontSize(12,0.01); hold on
% Output and compare data and fitted function coefs
text(n/10,3.6*m/4,sprintf('\t \t \t \t \t Amplitude \t X-Coord \t X-Width \t Y-Coord \t Y-Width \t Angle'),'Color','white');
text(n/10,3.8*m/4,sprintf('Fit \t %1.3f \t %1.0f \t %1.0f \t %1.0f \t %1.0f \t %1.2f',A),'Color','white');

%% Plot fitted data
subplot(3,8,[3 4 11 12]);
imagesc(x(1,:),y(:,1),fullGaussianFit); axis square;
colormap('gray');colorbar
title('Fitted gaussian')
upFontSize(12,0.01); hold on

%% Plot fitted data & ROI crop
maskedGaussianFit=fullGaussianFit.*mask;
subplot(3,8,[5 6 13 14]);
imagesc(x(1,:),y(:,1),maskedGaussianFit); axis square;
colormap('gray');colorbar;
title('Fitted gaussian with ROI overlay')
upFontSize(12,0.01); hold on

[centerX,centerY]=find(maskedGaussianFit==max(maskedGaussianFit(:)));
centerCoords=[centerX,centerY];

%% Plot thresholded contour plot
%nContourLevels=dataStruct.gaussianContourLevelMax(blockID);
subplot(3,8, [7 8 15 16]);
% plot underlying ROI masked gaussian response
fused=imfuse(maskedGaussianFit,bitmap,'ColorChannels',[2 1 1]);imagesc(x(1,:),y(:,1),fused);axis square; hold on;

% ! Changed from flipud on 8/21/2023 !
contourlinewidth=.5;
[mat,con]=contour((maskedGaussianFit),nContourLevels,'color',[1 1 1],'LineWidth',contourlinewidth,'ShowText','off'); %contour(flipud(maskedGaussianFit),nContourLevels,'red','LineWidth',.5,'ShowText','on');
axis square
contourLevels=con.LevelList;
upFontSize(12,0.01); hold on
title('Contour plot for generating masks')
clabel(mat,con,numel(con),'LabelSpacing',100)
con.LevelList = round(con.LevelList, 3);
for i=1:nContourLevels
    ROIMaskgaussian(i).threshPercentile=contourLevels(i);
    ROIMaskgaussian(i).area=maskedGaussianFit>ROIMaskgaussian(i).threshPercentile;
end
%colormap('gray');colorbar

%% Plot vertical and horizontal axes
vx_h=x(1,:); vy_v=y(:,1);
switch FitOrientation
    case 'fit', M=-tan(A(6));
        % generate points along _horizontal & _vertical axis
        vy_h=M*(vx_h-A(2))+A(4); hPoints = interp2(x,y,desiredResp,vx_h,vy_h,InterpMethod);
        vx_v=M*(A(4)-vy_v)+A(2); vPoints = interp2(x,y,desiredResp,vx_v,vy_v,InterpMethod);
    case 'dont', A(6)=0; 
        % generate points along _horizontal & _vertical axis
        vy_h=A(4)*ones(size(vx_h)); hPoints = interp2(x,y,desiredResp,vx_h,vy_h,InterpMethod);
        vx_v=A(2)*ones(size(vy_v)); vPoints = interp2(x,y,desiredResp,vx_v,vy_v,InterpMethod);
end
% plot lines 
hold on; plot(A(2),A(4),'+b',vx_h,vy_h,'r--',vx_v,vy_v,'b--','linewidth',2); box off; hold on
ylim([1 n]); xlim([1 m])

%% Plot long and short axes
dmin=1.1*min(desiredResp(:)); xposh = (vx_h-A(2))/cos(A(6))+A(2); xfit=xh(1,:); hfit=A(1)*exp(-(xfit-A(2)).^2/(2*A(3)^2));
dmax=1.1*max(desiredResp(:)); xposv = (vy_v-A(4))/cos(A(6))+A(4); yfit=yh(:,1); vfit=A(1)*exp(-(yfit-A(4)).^2/(2*A(5)^2));
subplot(3,8,[17:20]); 
plot(xposh,hPoints,'r',xfit,hfit,'k.'); grid on; box off;
title('Fitted and raw FFT amplitude along red axis')
upFontSize(12,0.005); hold on
subplot(3,8,[21:24]); 
plot(xposv,vPoints,'b',yfit,vfit,'k.'); grid on; box off;
title('Fitted and raw FFT amplitude along blue axis')
upFontSize(12,0.005); hold off
%set(gca,'YDir','reverse');
[~,h]=suplabel(['Fitted 2D gaussian and contour masks (response to small gaussian)'],'t',[.08 .08 .84 .87]);
if plotFlag==0
    close
end
end