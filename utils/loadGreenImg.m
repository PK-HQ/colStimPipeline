function [imgReference,imgTarget]=loadGreenImg(currentBlockStruct,referenceBlockStruct)
%% Find and load the reference and target green images
%{
filenameReferenceY=[mainPath  dataStructReference.monkey '/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameReferenceOI=['D:/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameTargetY=[mainPath  dataStructCurrent.monkey '/' dataStructCurrent.monkey dataStructCurrent.date '/green' num2str(dataStructCurrent.greenImgID) '_binned.bmp'];
filenameTargetOI=['D:/'  dataStructCurrent.monkey dataStructCurrent.date '/green' num2str(dataStructCurrent.greenImgID) '_binned.bmp'];
%}
imgReference=nan(512,512);
imgTarget=nan(512,512);
if isfile(referenceBlockStruct.greenServer)
    imgReference=imread(referenceBlockStruct.greenServer);
    %disp(referenceBlockStruct.greenServer)
else
    imgReference=imread(referenceBlockStruct.greenOI);
    %disp(referenceBlockStruct.greenOI)
end
if isfile(currentBlockStruct.greenServer)
    imgTarget=imread(currentBlockStruct.greenServer);
    %disp(currentBlockStruct.greenServer)
elseif isfile(currentBlockStruct.greenOI)
    imgTarget=imread(currentBlockStruct.greenOI);
    %disp(currentBlockStruct.greenOI)
end
%{
imgTarget=imread(imgTargetName);
disp(imgTargetName)
imgReference=imread(imgReferenceName);
disp(imgReferenceName)
%}
end