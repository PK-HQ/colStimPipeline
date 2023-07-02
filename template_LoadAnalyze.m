function template_LoadAnalyze(dataStruct, plotOrNot)

%% What it does
% Loads file specified in dataStruct.field 
% Run whatever you want
%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Get filenames
filenameStruct=generateFilenames(dataStruct); %generates filenames for monkey/run/session in dataStruct
%pdfFilename=filenameStruct.colmapPDF; % here is the save filename

%% Load selected variables from filename
load(filenameStruct.Orientation,'Ort','RespCondPCA','Mask','MapOrt','MapAmpOrt','ColorMap'); %contains mask, RespCondPCA
%conds
end