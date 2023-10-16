function index=findMatch(struct,regexString)
logicalIdx = ~cellfun(@isempty,regexp(struct,regexString));
index=find(logicalIdx==1);
end