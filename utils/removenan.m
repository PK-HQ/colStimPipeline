function dataCleaned=removenan(data)
data_flat = data(:);

% Find the indices of non-NaN elements
nonNanIndices = ~isnan(data_flat);

% Extract the non-NaN elements using logical indexing
dataCleaned = data_flat(nonNanIndices);
end