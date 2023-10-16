function firstfilename=getFirstFile(folderpath,regex)
filenames=dir([folderpath regex]);
if ~isempty(filenames)
    firstfilename=[folderpath filenames(1).name];
end
end