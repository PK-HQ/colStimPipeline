function dataMat=struct2mat(dataStruct)
    dataMat=cell2mat(struct2cell(dataStruct));
end