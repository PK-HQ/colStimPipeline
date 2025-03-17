function B = FilterFermi3D(image, LowCutOff, HighCutOff, SizePxl)

% Check input/output arguments
SizeImage = size(image);
Height = SizeImage(1);
Width = SizeImage(2);
NumSlices = SizeImage(3); % Number of 2D slices

% Fermi filter parameters
ParmFermiLowPass = [1, HighCutOff, 0, HighCutOff * 0.05];
ParmFermiHighPass = [1, LowCutOff, 0, LowCutOff * 0.05];

% Generate spatial frequency grid
SFX = ((1:Width) - floor(Width / 2) - 1) / (Width - 1) / SizePxl;
SFY = ((1:Height) - floor(Height / 2) - 1) / (Height - 1) / SizePxl;
[SFXX, SFYY] = meshgrid(SFX, SFY);
SF2D = abs(SFXX + 1i * SFYY);

% Create Fermi filters
if HighCutOff == inf
    FiltFermiLowPass = zeros(Height, Width);
else
    FiltFermiLowPass = FuncWoNFermi(ParmFermiLowPass, SF2D);
end
if LowCutOff == 0
    FiltFermiHighPass = ones(Height, Width);
else
    FiltFermiHighPass = FuncWoNFermi(ParmFermiHighPass, SF2D);
end
FiltFermi = FiltFermiHighPass - FiltFermiLowPass;

% Expand filter to 3D for processing all slices
FiltFermi3D = repmat(FiltFermi, 1, 1, NumSlices);

% Apply the filter to the 3D image
B = ifft(ifft(ifftshift(ifftshift( ...
    fftshift(fftshift(fft(fft(image, [], 1), [], 2), 1), 2) .* ...
    FiltFermi3D, 1), 2), [], 1), [], 2);

% Ensure real output if the input is real
if isreal(image)
    B = real(B);
end
end
