%% Fit a 2D gaussian function to data
%% PURPOSE:  Fit a 2D gaussian centroid to simulated data
% Uses lsqcurvefit to fit
%
% INPUT:
% 
%   MdataSize: Size of nxn data matrix
%   x0 = [Amp,x0,wx,y0,wy,fi]: Inital guess parameters
%   x = [Amp,x0,wx,y0,wy,fi]: simulated centroid parameters
%   noise: noise in % of centroid peak value (x(1)
%   InterpolationMethod = 'nearest' or 'linear' or 'spline' or 'cubic'
%       used to calculate cross-section along minor and major axis
%     
%
%
% NOTE:
%   The initial values in x0 must be close to x in order for the fit
%   to converge to the values of x (especially if noise is added)
%
% OUTPUT:  non
%
% CREATED: G. Nootz  May 2012
% 
%  Modifications:
%  non

function fit2Dgaussian(varargin)
if length(varargin)==1
    desiredResp=varargin{1};
    mask=ones(size(desiredResp));
elseif length(varargin)==2
    desiredResp=varargin{1};
    mask=varargin{2};
end
%% ---------User Input---------------------
MdataSize = size(desiredResp,1); % Size of nxn data matrix
% parameters are: [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
x0 = [1,200,100,250,100,0.7854]; %Inital guess parameters
noise = 10; % noise in % of centroid peak value (x(1))
InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
FitForOrientation = 1; % 1: fit for orientation. 0: do not fit for orientation

%% ---Generate centroid to be fitted--------------------------------------
[X,Y] = meshgrid(1:MdataSize);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
[Xhr,Yhr] = meshgrid(linspace(1,512,300)); % generate high res grid for plot
xdatahr = zeros(300,300,2);
xdatahr(:,:,1) = Xhr;
xdatahr(:,:,2) = Yhr;
Z = desiredResp;


%% Fit; x = [Amp, x0, wx, y0, wy, fi]
if FitForOrientation == 0
    x0 =x0(1:5);
    xin(6) = 0; 
    x =x(1:5);
    lb = [0,1,0,1,0];
    ub = [realmax('double'),MdataSize,MdataSize^2,MdataSize,MdataSize^2];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub);
    x(6) = 0;
else
    % define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
    lb = [0,-MdataSize/2,0,-MdataSize/2,0,-pi/4];
    ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2,pi/4];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z,lb,ub);
end


% Fitted params
headerStr = ['             Amplitude','    x_0', '    \sigma_x','    y_0','    \sigma_y','     \theta'];
fittedParamStr = ['Fitted      ',num2str(x(1), '% 100.4f'),'       ',num2str(x(2), '% 100.0f'),'   ',num2str(x(3), '% 100.0f'),'   ',num2str(x(4), '% 100.0f'),'   ',num2str(x(5), '% 100.0f'),'   ',num2str(x(6), '% 100.2f')];

%% ---------Plot 3D Image-------------


%% -----Plot profiles----------------
hf2 = figure(1);
subplot(4,4, [4])
C = del2(Z);
mesh(X,Y,Z,C) %plot data
hold on
surface(Xhr,Yhr,D2GaussFunctionRot(x,xdatahr),'EdgeColor','none') %plot fit
%axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5 -noise noise+x(1)])
alpha(0.2)  
hold on

set(hf2, 'Position', [20 20 950 900])
alpha(0)
subplot(4,4, [5,6,7,9,10,11,13,14,15])
imagesc(X(1,:),Y(:,1)',Z)
set(gca,'YDir','reverse')
colormap('gray')
% Text for printed params
text(MdataSize/2*0,MdataSize/4*3.2,headerStr,'Color', 'yellow','Interpreter', 'tex')
text(MdataSize/2*0,MdataSize/4*3.4,fittedParamStr,'Color','yellow')





g1d7*g1d7'


%% -----Calculate cross sections-------------
% generate points along horizontal axis
m = -tan(x(6));% Point slope formula
b = (-m*x(2) + x(4));
xvh = 1:512; %-MdataSize/2:MdataSize/2;
yvh = xvh*m + b;
hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
% generate points along vertical axis
mrot = -m;
brot = (mrot*x(4) - x(2));
yvv = 1:512;%-MdataSize/2:MdataSize/2;
xvv = yvv*mrot - brot;
vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);

hold on % Indicate major and minor axis on plot

% % plot pints 
% plot(xvh,yvh,'r.') 
% plot(xvv,yvv,'g.')

% plot lins 
plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r');box off 
plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'b');box off
upFontSize(13,.015)

hold off
%axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5])

%% Plot 1D gaussian profiles
ymin = - noise * x(1);
ymax = x(1)*(1+noise);
xdatafit = linspace(1,512,300);
hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
vdatafit = x(1)*exp(-(xdatafit-x(4)).^2/(2*x(5)^2));
subplot(4,4, [1:3])
xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black');box off;upFontSize(13,.015)

subplot(4,4,[8,12,16])
xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
plot(vPoints,xposv,'b.',vdatafit,xdatafit,'black');box off;upFontSize(13,.015)
set(gca,'YDir','reverse')
figure(gcf) % bring current figure to front
1;

upFontSize(13,.015)
suplabel('2D gaussian fit','t')

end



 
