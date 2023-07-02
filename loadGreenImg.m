function [imgReference,imgTarget]=loadGreenImg(dataStructReference,dataStructCurrent)
%% Find and load the reference and target green images
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end

filenameReferenceY=[mainPath  dataStructReference.monkey '/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameReferenceOI=['D:/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameTargetY=[mainPath  dataStructCurrent.monkey '/' dataStructCurrent.monkey dataStructCurrent.date '/green2_binned.bmp'];
filenameTargetOI=['D:/'  dataStructCurrent.monkey dataStructCurrent.date '/green2_binned.bmp'];
if isfile(filenameReferenceY)
    imgReference=imread(filenameReferenceY);
else
    imgReference=imread(filenameReferenceOI);
end
if isfile(filenameTargetY)
    imgTarget=imread(filenameTargetY);
elseif isfile(filenameTargetOI)
    imgTarget=imread(filenameTargetOI);
end
end