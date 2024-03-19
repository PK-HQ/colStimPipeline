function [Mask, RespCondPCA, Ort]=loadOrtMapData(referenceEntry)
    Mask=[];
    RespCondPCA=[]; 
    Ort=[];
    if isfile(referenceEntry.Orientation)
        load(referenceEntry.Orientation, 'Mask', 'RespCondPCA', 'Ort');
    end
end