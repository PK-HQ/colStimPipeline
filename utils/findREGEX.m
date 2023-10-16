function [fullFileNames]=findREGEX(regexStr)
  fullFileNames={};
  % get struct of all matches within current directory
  filesStruct=dir(regexStr);
  % if not empty, extract filename and folder path
  if ~isempty(filesStruct)
    folderPath=extractfield(filesStruct,"folder");
    fileNames=extractfield(filesStruct,"name");
    % make full filenames
    for fileNo=1:length(fileNames)
      fullFileNames{fileNo}=fullfile(folderPath{fileNo},fileNames{fileNo});
    end
  end
end