function savePDF(filename, monkeyName, saveFlag, plotNo)
%% Saves plot as PDF
% filename is a string
% monkeyname is a string for monkey name
% saveFlag == 1 means save
% plotNo == 1 means save as a new file, >1 means append
filenameMonkey=['Y:/' monkeyName '/Meta/' filename '.pdf'];
filenameMeetings=['Y:/users/PK/Eyal/meetings/' filename '.pdf'];
if plotNo==1 && saveFlag==1
    export_fig(filenameMonkey,'-pdf','-nocrop');
    copyfile(filenameMonkey, filenameMeetings)
elseif plotNo>1 && saveFlag==1
    export_fig(filenameMonkey,'-pdf','-append','-nocrop');
    copyfile(filenameMonkey, filenameMeetings)
end
end