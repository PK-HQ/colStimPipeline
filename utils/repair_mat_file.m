function fixedData = repair_mat_file(inputFilename)
%REPAIR_MAT_FILE Repair a MAT v5 file with a bad subsystem data pointer.
%
% fixedData = repair_mat_file(inputFilename)
%
% This creates inputName_fixed.mat in the same folder as inputFilename,
% clears the subsystem-data offset in the MAT-file header, loads the fixed
% file, and returns the loaded struct. inputFilename is never modified.

if nargin ~= 1 || isempty(inputFilename)
    error('Provide a path to the damaged .mat file.');
end

inputFilename = char(inputFilename);
[folder, base, ext] = fileparts(inputFilename);
if isempty(ext)
    ext = '.mat';
end
fixedFilename = fullfile(folder, [base '_fixed' ext]);

if exist(fixedFilename, 'file')
    error('Fixed file already exists, so it was not overwritten: %s', fixedFilename);
end

info = dir(inputFilename);
if isempty(info)
    error('Input file does not exist: %s', inputFilename);
end
fileSize = info.bytes;

fid = fopen(inputFilename, 'rb');
if fid < 0
    error('Could not open input file: %s', inputFilename);
end
cleanupIn = onCleanup(@() fclose(fid));

header = fread(fid, 128, 'uint8=>uint8')';
if numel(header) ~= 128
    error('File is too small to be a MAT v5 file: %s', inputFilename);
end

headerText = char(header(1:19));
if ~strcmp(headerText, 'MATLAB 5.0 MAT-file')
    error('This repair is only intended for MATLAB 5.0 MAT-files.');
end

endianKey = char(header(127:128));
if ~strcmp(endianKey, 'IM')
    error('This repair currently supports little-endian MAT files with endian key IM.');
end

subsystemOffset = double(typecast(uint8(header(117:124)), 'uint64'));
if subsystemOffset <= 128 || subsystemOffset > fileSize
    error('No valid trailing subsystem stream was found. subsystemOffset=%g, fileSize=%g.', subsystemOffset, fileSize);
end

fprintf('Input:  %s\n', inputFilename);
fprintf('Output: %s\n', fixedFilename);
fprintf('Clearing subsystem offset %.0f and preserving %.0f file bytes.\n', ...
    subsystemOffset, fileSize);

header(117:124) = 0;

outFid = fopen(fixedFilename, 'wb');
if outFid < 0
    error('Could not create output file: %s', fixedFilename);
end
cleanupOut = onCleanup(@() fclose(outFid));

fwrite(outFid, header, 'uint8');
fseek(fid, 128, 'bof');

remaining = fileSize - 128;
bufferSize = 1024 * 1024;
while remaining > 0
    count = min(bufferSize, remaining);
    chunk = fread(fid, count, 'uint8=>uint8');
    if numel(chunk) ~= count
        error('Unexpected end of file while copying valid MAT stream.');
    end
    fwrite(outFid, chunk, 'uint8');
    remaining = remaining - count;
end

clear cleanupOut cleanupIn

fprintf('\nChecking repaired file with whos -file ...\n');
vars = whos('-file', fixedFilename);
disp({vars.name}');

fixedData = load(fixedFilename);
fprintf('Repair complete.\n');
end
