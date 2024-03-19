function [Mask, RespCondPCA, Ort] = loadReferenceData(filenameStructReference)
    % Load mask and response condition data
    load(filenameStructReference.Orientation, 'Mask', 'RespCondPCA', 'Ort');
end