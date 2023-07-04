function extractPCA(selectedOrt,RespCondPCA,ROImaskNaN,ROISquareMask)
responsePCA=RespCondPCA(:,:,selectedOrt);
for ortNo=1:size(responsePCA,3)
    %VER for orientation, masked with SNR mask
    VERpca(:,:,ortNo)=responsePCA(:,:,ortNo).*ROImaskNaN;

    %VER for orientation, masked with ROI mask (to isolate site)
    VERpca(:,:,ortNo)=VERpca(:,:,ortNo).*ROISquareMask;
end
end