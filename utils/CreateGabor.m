%% Create Gabor BMP files
% 1. Spatial frequency
% 2. Bandwidth
%
% Abbreviation in file names
% Z---size
% S---spatial frequency
% W--bandwidth
% O---orientation
% P---phase
%
%
% YC at ES lab
% Created on Jun. 10, 2004
% Last modified on May 12, 2008

%% Clear workspace
clear all;
close all;

%% Select a directory
PathName = uigetdir('.','Select a Directory to Save BMP Files');
drawnow;

if PathName==false
  disp('No directory is selected!');
  return;
end

%% Parameters
RatioPxlPerDeg = 50;  % pixels/deg, ratio of screen to visul field
GaborSizes=[1.6:.1:1.9];
GaborOrts=[0 90];%[0:10:180,22.5,45,67.5,112.5,135,157.5];
%% Create BMP
for GaborSize = GaborSizes % deg, size (square)
  % Some parameters and arrays
  GaborSDX = GaborSize*RatioPxlPerDeg/6;  % pixels, sigma X, 1/6 of whole size
  GaborSDY = GaborSDX;  % pixels, sigma Y

  SizeX = ceil(GaborSize*RatioPxlPerDeg/2)*2+1;
  SizeY = SizeX;

  X = -(SizeX-1)/2:(SizeX-1)/2;
  Y = -(SizeY-1)/2:(SizeY-1)/2;

  % Blank
  A = ones(length(Y),length(X))*0.5;
  FileName = ...
    sprintf('BlankZ%04g.bmp',GaborSize*100);
  imwrite(A,fullfile(PathName,FileName),'bmp');

  % Gabor
  for GaborOrt = GaborOrts  % deg, orientation
    for GaborPhs = 0:90:270  % deg, phase
      % Spatial frequency
      for GaborSF = 0.25*2.^(0:7)  % oct, band width
        A = Gabor2D(X,Y,[GaborOrt,GaborSF/RatioPxlPerDeg, ...
                         GaborSDX,GaborSDY,GaborPhs]);
        A = A/max(abs(A(:)))*0.5+0.5;  % nomalize A to [0,1]

        FileName = ...
          sprintf('GaborZ%04gS%04dO%05dP%05d.bmp', ...
                  GaborSize*100,GaborSF*100,GaborOrt*100,GaborPhs*100);
        imwrite(A,fullfile(PathName,FileName),'bmp');
      end
      % Bandwidth
      for GaborBW = 0.25:0.25:2  % oct, band width
        GaborSF = ...
          1/pi*sqrt(log(2)/2)*(2^GaborBW+1)/(2^GaborBW-1)/ ...
          GaborSDX*RatioPxlPerDeg;

        A = Gabor2D(X,Y,[GaborOrt,GaborSF/RatioPxlPerDeg, ...
                         GaborSDX,GaborSDY,GaborPhs]);
        A = A/max(abs(A(:)))*0.5+0.5;  % nomalize A to [0,1]

        FileName = ...
          sprintf('GaborZ%04gW%04dO%05dP%05d.bmp', ...
                  GaborSize*100,GaborBW*100,GaborOrt*100,GaborPhs*100);
        imwrite(A,fullfile(PathName,FileName),'bmp');
      end
    end
  end
end


