function argOut = fastload(pathToMatfile, myVarName)
% Get the location of the variable in the file using hdf syntax
% / by itself is the root of the file, then variables names come after
% Note: also works very nicely for nested structures where a.b.c.d.e would
% have varName as /a/b/c/d/e. 
h5loc = ['/' myVarName]; % Always /, not like windows/linux filesep
% Open the file using H5F. 
fid = H5F.open(pathToMatfile); 
% Open the file using H5D. 
dsetid = H5D.open(fid,h5loc); 
% Load in the dataset
argOut = H5D.read(dsetid,h5loc); % All done
% Clean up
H5D.close(dsetid);
H5F.close(fid);
end