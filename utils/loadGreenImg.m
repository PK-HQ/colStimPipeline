function [imgReference,imgTarget]=loadGreenImg(dataStructReference,dataStructCurrent,imgTargetName,imgReferenceName)
%% Find and load the reference and target green images
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end

filenameReferenceY=[mainPath  dataStructReference.monkey '/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameReferenceOI=['D:/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameTargetY=[mainPath  dataStructCurrent.monkey '/' dataStructCurrent.monkey dataStructCurrent.date '/green' num2str(dataStructCurrent.greenImgID) '_binned.bmp'];
filenameTargetOI=['D:/'  dataStructCurrent.monkey dataStructCurrent.date '/green' num2str(dataStructCurrent.greenImgID) '_binned.bmp'];
if isfile(filenameReferenceY)
    imgReference=imread(filenameReferenceY);
    disp(filenameReferenceY)
else
    imgReference=imread(filenameReferenceOI);
    disp(filenameReferenceOI)
end
if isfile(filenameTargetY)
    imgTarget=imread(filenameTargetY);
    disp(filenameTargetY)
elseif isfile(filenameTargetOI)
    imgTarget=imread(filenameTargetOI);
    disp(filenameTargetOI)
end
%{
imgTarget=imread(imgTargetName);
disp(imgTargetName)
imgReference=imread(imgReferenceName);
disp(imgReferenceName)
%}
end