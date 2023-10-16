function bitmap=createDCbitmap(columnsThinned,DC)

if DC==0.1
    dimDCpix=[15,15]; %shoudl calculate from projector DMD dims
end

pixels=numel(find(columnsThinned>0));
[xPositions,yPositions]=find(columnsThinned>0);
for pixel=1:pixels
    xPos=xPositions(pixel);
    yPos=yPositions(pixel);
    minXBound=(dimDCpix(1)-1)/2;
    maxXBound=(dimDCpix(1)+1)/2;
    minYBound=(dimDCpix(2)-1)/2;
    maxYBound=(dimDCpix(2)+1)/2;
    % place a DC=0.1 square centered on the pixel
    columnsThinned(xPos-minXBound:xPos+maxXBound,...
        yPos-minYBound:yPos+maxYBound)=1; %check if this is round number
    bitmap=columnsThinned>0;
end
end