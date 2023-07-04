function normVERpca=normalizePCA(VERpca,ROImaskNaN)
for ortNo=1:size(VERpca,3)
  map=VERpca(:,:,ortNo);
  normMap=rescale(map,0,1); 
  %clips image, remap to -1 to 1
  normMap = normMap .* ROImaskNaN;
  normVERpca(:,:,ortNo)=normMap;
end
end