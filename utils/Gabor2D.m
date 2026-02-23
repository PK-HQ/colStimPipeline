function Z = Gabor2D(X,Y,Parameters)

% Two dimensional Gabor function
% Parameters = [GaborOrt,GaborSF,GaborSDX,GaborSDY,GaborPhs];

GaborOrt = Parameters(1)/180*pi;
GaborSF  = 2*pi*Parameters(2);
GaborSDX = Parameters(3);
GaborSDY = Parameters(4);
GaborPhs = Parameters(5)/180*pi;

[XX,YY] = meshgrid(X,Y);

Z = cos(GaborSF*sin(GaborOrt)*XX+GaborSF*cos(GaborOrt)*YY+GaborPhs).* ...
    exp(-(XX.*XX)/(2*GaborSDX*GaborSDX)-(YY.*YY)/(2*GaborSDY*GaborSDY))/ ...
    (2*pi*GaborSDX*GaborSDY);


