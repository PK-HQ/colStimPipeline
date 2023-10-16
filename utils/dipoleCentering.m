TS.Header.Conditions.PosCond(1:2)

load('Y:\Chip\Chip20220909\dipole\run3\M28D20220909R2StabDFFTAmpS004E023PF0400.mat')
-1.5000   -0.5500
load('Y:\Chip\Chip20220909\dipole\run4\M28D20220909R3StabDFFTAmpS004E023PF0400.mat')
-2.2000   -0.5000
load('Y:\Chip\Chip20220909\dipole\run5\M28D20220909R4StabDFFTAmpS004E023PF0400.mat')
-2.3500   -0.7000

DataCondBlank=DataCond(:,:,1);
DataCondSignal=DataCond(:,:,2:end) - DataCondBlank;
figure; colormap(fireice)
subplot(2,1,1)
imgsc(DataCondSignal(:,:,1)-DataCondSignal(:,:,2))
subplot(2,1,2)
imgsc(DataCondSignal(:,:,3)-DataCondSignal(:,:,4))


subplot(2,1,2)
imgsc(DataCondSignal(:,:,2)-DataCondSignal(:,:,1))