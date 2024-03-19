function [filenameStructCurrent,filenameStructReference]=grabFilenames(dsCurrentSess, dsReferenceSess)
    filenameStructCurrent = generateFilenames(dsCurrentSess);
    filenameStructReference = generateFilenames(dsReferenceSess);
end