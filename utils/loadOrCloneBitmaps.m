function [bmp,status] = loadOrCloneBitmaps(currentBlockStruct, verbose)
% LOADORCLONEBITMAPS  Load stimulus bitmaps for a currentBlockStruct; clone if necessary.
%
%   [bmp,status] = loadOrCloneBitmaps(currentBlockStruct)
%
% Inputs
%   currentBlockStruct   : struct with fields
%                .TS       (full path to .mat containing TS variable)
%                .ROITC    (destination folder for bitmaps)
%                .mainPath (root repo that holds users/PK/colStimPipeline)
%                .date     (YYYY-MM-DD or similar sub-folder under V:/PK/ColSeries)
%   verbose : optional logical (default=true).  Print progress messages.
%
% Outputs
%   bmp     : 1080×1920×N double.  NaN sheets where file missing.
%   status  : struct with fields .loaded, .cloned, .missing (logical vectors)
%
% Requires   findMatchingFiles.m  (same helper you already have)
%
% -------------------------------------------------------------------------

if nargin<2, verbose = true; end
bmp = []; status = struct('loaded',[],'cloned',[],'missing',[]);

%% ---------- 1. Sanity checks -------------------------------------------------
if ~isfile(currentBlockStruct.TS)
    warning('TS file %s not found.', currentBlockStruct.TS);
    bmp = nan(1080,1920,2);    % keep interface identical
    status.missing = true(1,2);
    return
end

%% ---------- 2. Parse TS for bitmap names ------------------------------------
TSstruct  = load(currentBlockStruct.TS,"TS");
TS=TSstruct.TS;
proj = TS.Header.Conditions.ProjImg;
proj = proj(cellfun('isempty',strfind(proj,'Dot')));  % drop “Dot” stimuli
proj = unique(proj);
if numel(proj)>2                      % keep last two by convention
    proj = proj(end-1:end);
end
proj = strrep(proj,'\','/');          % normalise slashes

% 1. Call fileparts and get all names and extensions separately
[~, names, exts] = cellfun(@fileparts, proj, 'UniformOutput', false); 
% 2. Concatenate the names and extensions
rawNames = cellfun(@strcat, names, exts, 'UniformOutput', false);
nImgs    = numel(rawNames);

%% ---------- 3. Attempt 1 – LOAD from ROITC ----------------------------------
[bmp, found] = tryLoad(currentBlockStruct.ROITC, rawNames, nImgs);
if all(found)
    status.loaded  = found;
    status.cloned  = false(1,nImgs);
    status.missing = false(1,nImgs);
    if verbose, fprintf('Loaded %d/%d bitmaps from ROITC.\n',nImgs,nImgs); end
    return
end

%% ---------- 4. Attempt 2 – CLONE missing files then LOAD again --------------
if verbose
    need = find(~found);
    fprintf('Missing %d bitmap(s); cloning those now...\n',numel(need));
end

% Build array of repo-relative names (after "PK")
repoRel = cellfun(@(p) extractAfter(p,'PK'), proj,'Uni',false);

clonePaths = findMatchingFiles( ...
              fullfile(currentBlockStruct.mainPath,'users','PK','colStimPipeline'), ...
              repoRel, currentBlockStruct.date);

% ensure destination exists
if ~exist(currentBlockStruct.ROITC,'dir'), mkdir(currentBlockStruct.ROITC); end

for i = 1:numel(clonePaths)
    src = clonePaths{i};
    if isempty(src) || ~isfile(src); continue; end   % skip if search failed
    [~,nm,ext] = fileparts(src);
    dst = fullfile(currentBlockStruct.ROITC,[nm ext]);
    try
        copyfile(src,dst); if verbose, fprintf('Copied %s\n',[nm ext]); end
    catch ME
        warning('Could not copy %s ? %s (%s)',src,dst,ME.message);
    end
end

% second load attempt
[bmp2,found2] = tryLoad(currentBlockStruct.ROITC, rawNames, nImgs);

% merge results: keep already-loaded frames, update cloned ones
bmp(:,:,~found) = bmp2(:,:,~found);

status.loaded  = found|found2;
status.cloned  = ~found & found2;
status.missing = ~status.loaded;

if verbose
    fprintf('Final: %d loaded, %d cloned, %d missing.\n', ...
            nnz(status.loaded & ~status.cloned), ...
            nnz(status.cloned), nnz(status.missing));
end

end  % ----------- main --------------------------------------------------------

% =====================================================================
function [imgStack, ok] = tryLoad(folder, names, nImgs)
% Load bitmaps 'names' from 'folder' ? imgStack (1080×1920×N)
imgStack = nan(1080,1920,nImgs);
ok       = false(1,nImgs);

for k = 1:nImgs
    f = fullfile(folder, names{k});
    if exist(f,'file')
        imgStack(:,:,k) = double(imread(f));
        ok(k) = true;
    end
end
end
