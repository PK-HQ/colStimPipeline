function createSuperblock(filenamestruct)
superblockFlag=~isempty(datastruct(entryNo).baselinerun);
switch superblockFlag
    case {1}
        %load TS.mat and FFT.mat for baseline
        
        %load TS.mat and FFT.mat for columnar optostim runs
        
        % merge TS file into the format of the regular combined blocks (OR identify conditions based on the condtable)
        
        % merge FFT.mat file
        
        %output
        
        
    case {0}
        % do nothing
end