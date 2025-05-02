load('E:\Pepper20250328\run1\Data_vdaqlog.mat')
respFFT=real(FFTCond);
respFFTdipole=respFFT(:,:,[2 3])-respFFT(:,:,1);
respFFTdipole0=respFFTdipole(:,:,1);
respFFTdipole90=respFFTdipole(:,:,2);
subplot(1,3,1); title('0 dipole')
imgsc(respFFTdipole0); colorbar
axis square; title('0 dipole')
subplot(1,3,2); 
imgsc(respFFTdipole90);colorbar
axis square; title('90 dipole')
subplot(1,3,3);
imgsc(respFFTdipole90-respFFTdipole0);colorbar
axis square; title('90-0 dipole')
