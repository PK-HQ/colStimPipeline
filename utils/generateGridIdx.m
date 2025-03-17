function plotIdx = generateGridIdx(row,column)
plotIdx=reshape(1:row*column,column,row)';
end
