for run=6:11%1:6
  subplot(3,5,run)
  
  load(['Y:\Chip\Chip20230124\run' num2str(run) '\M28D20230124R' num2str(run) 'StabDFFTAmpS004E023PF0400.mat'])
  load(['Y:\Chip\Chip20230124\run' num2str(run) '\M28D20230124R' num2str(run) 'TS.mat'])
  blank=DataCond(:,:,1);
  signal=DataCond(:,:,2:3)-blank;

  
  diffSignal=signal(:,:,2) - signal(:,:,1);
  imgsc(imgaussfilt(diffSignal,2))
  title(TS.Header.Conditions.PosCond(1:2))
  colormap(gray)
end