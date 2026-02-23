function getBitmapPipelineCode(zipPath, mode, outFolder, varargin)
% package_mainPipeline_mode  Build a minimal code bundle for mainPipeline.m for a given analysisMode.
%
% It:
%   1) Unzips zipPath into a temp folder
%   2) Parses mainPipeline.m to find:
%        - preamble (before "switch analysisMode")
%        - the "case '<mode>'" block
%   3) Collects functions/scripts referenced in those sections that exist inside the project
%      (so MATLAB built-ins/toolbox functions are excluded)
%   4) Uses matlab.codetools.requiredFilesAndProducts to recursively get dependencies
%   5) Copies only those files into outFolder, preserving folder structure
%
% Usage:
%   package_mainPipeline_mode('colStimPipeline.zip', 'expt', 'bitmapPipeline')
%   package_mainPipeline_mode('colStimPipeline.zip', 'summary', 'bundle_summary')
%
% Optional name-value args:
%   'ExtraEntryPoints' : cellstr of additional .m files or function names to seed dependency search
%   'IncludeScriptsWithPattern' : cellstr patterns to include (e.g., {'exptListBiasingFull'})
%   'Verbose' : true/false
%
% Notes:
% - This is static analysis + dependency closure of referenced project functions.
% - If your code uses eval/feval/function handles built from strings, you may need
%   to add those as ExtraEntryPoints or pattern includes.

p = inputParser;
p.addRequired('zipPath', @(x)ischar(x)||isstring(x));
p.addRequired('mode', @(x)ischar(x)||isstring(x));
p.addRequired('outFolder', @(x)ischar(x)||isstring(x));
p.addParameter('ExtraEntryPoints', {}, @(x)iscell(x) || isstring(x));
p.addParameter('IncludeScriptsWithPattern', {'exptListBiasingFull'}, @(x)iscell(x) || isstring(x));
p.addParameter('Verbose', true, @(x)islogical(x) && isscalar(x));
p.parse(zipPath, mode, outFolder, varargin{:});

zipPath  = char(p.Results.zipPath);
mode     = char(p.Results.mode);
outFolder = char(p.Results.outFolder);
extraEntry = cellstr(p.Results.ExtraEntryPoints);
patternIncludes = cellstr(p.Results.IncludeScriptsWithPattern);
verbose = p.Results.Verbose;

if ~exist(zipPath, 'file')
    error('Zip not found: %s', zipPath);
end

tmpDir = tempname;
mkdir(tmpDir);
unzip(zipPath, tmpDir);

mainFile = fullfile(tmpDir, 'mainPipeline.m');
if ~exist(mainFile, 'file')
    error('mainPipeline.m not found at zip root. Found files under: %s', tmpDir);
end

% Index all .m files in project (relative path)
mFiles = dir(fullfile(tmpDir, '**', '*.m'));
relPaths = cell(numel(mFiles),1);
baseNames = cell(numel(mFiles),1);
absPaths = cell(numel(mFiles),1);
for i = 1:numel(mFiles)
    absPaths{i} = fullfile(mFiles(i).folder, mFiles(i).name);
    relPaths{i} = erase(absPaths{i}, [tmpDir filesep]);
    baseNames{i} = mFiles(i).name; % includes ".m"
end

% Map function/script base name -> abs path (handle collisions)
fileMap = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(mFiles)
    key = baseNames{i};
    if ~isKey(fileMap, key)
        fileMap(key) = absPaths{i};
    else
        % Collision: keep the shorter rel path (usually more "primary"), but warn.
        oldAbs = fileMap(key);
        oldRel = erase(oldAbs, [tmpDir filesep]);
        newRel = relPaths{i};
        if numel(newRel) < numel(oldRel)
            fileMap(key) = absPaths{i};
        end
    end
end

% Read mainPipeline.m
txt = fileread(mainFile);
lines = splitlines(txt);

% Find "switch analysisMode"
switchLine = find(contains(lines, 'switch analysisMode'), 1, 'first');
if isempty(switchLine)
    error('Could not find "switch analysisMode" in mainPipeline.m');
end

% Extract preamble text (before switch)
preambleText = strjoin(lines(1:switchLine-1), newline);

% Extract the requested case block
caseIdx = find(contains(lines, "case '" + string(mode) + "'"), 1, 'first');
if isempty(caseIdx)
    error("Mode '%s' not found in mainPipeline.m (no case '%s')", mode, mode);
end

% Determine where the case block ends: next line that starts with "case " or "otherwise"
endIdx = numel(lines);
for i = caseIdx+1:numel(lines)
    L = strtrim(lines(i));
    if startsWith(L, "case '") || startsWith(L, 'otherwise')
        endIdx = i-1;
        break;
    end
end
caseText = strjoin(lines(caseIdx:endIdx), newline);

% Collect candidate function calls from preamble+caseText.
scanText = preambleText + newline + caseText;

% Regex: functionName(
tokens = regexp(scanText, '\<([A-Za-z]\w*)\s*\(', 'tokens');
funcNames = unique(string([tokens{:}]));

% Drop obvious MATLAB keywords (not exhaustive, but fine)
keywords = ["if","for","while","switch","case","otherwise","end","return","function", ...
            "zeros","ones","nan","size","length","error","warning","disp","fprintf", ...
            "fullfile","dir","exist","isfield","ischar","isstring","strcmp","strcmpi", ...
            "cell","struct","double","single","logical","isempty","find","contains", ...
            "cat","reshape","permute","squeeze","mean","sum","std","min","max","round", ...
            "floor","ceil","interp1","interp2","imread","imwrite","imagesc","plot", ...
            "hold","title","xlabel","ylabel","legend","axis","set","get","gca","gcf"];

funcNames = funcNames(~ismember(funcNames, keywords));

% Keep only functions/scripts that exist as .m files in the project
entryFiles = strings(0,1);

% Always include mainPipeline.m itself
entryFiles(end+1) = mainFile;

for f = funcNames(:)'
    key = char(f + ".m");
    if isKey(fileMap, key)
        entryFiles(end+1) = string(fileMap(key));
    end
end

% Include script files matching patterns (helps with setupEnv-run scripts, etc.)
for pat = string(patternIncludes(:))'
    hits = endsWith(string(baseNames), ".m") & contains(string(baseNames), pat);
    entryFiles = [entryFiles; string(absPaths(hits))]; %#ok<AGROW>
end

% Include any extra entry points requested
for e = string(extraEntry(:))'
    if endsWith(e, ".m")
        % treat as relative path inside zip OR absolute path
        candidate = char(e);
        if exist(candidate, 'file')
            entryFiles(end+1) = string(candidate);
        else
            candidate2 = fullfile(tmpDir, candidate);
            if exist(candidate2, 'file')
                entryFiles(end+1) = string(candidate2);
            else
                % maybe user passed just "foo" without .m
                key = char(e + ".m");
                if isKey(fileMap, key)
                    entryFiles(end+1) = string(fileMap(key));
                else
                    warning('ExtraEntryPoints item not found: %s', e);
                end
            end
        end
    else
        key = char(e + ".m");
        if isKey(fileMap, key)
            entryFiles(end+1) = string(fileMap(key));
        else
            warning('ExtraEntryPoints item not found: %s', e);
        end
    end
end

entryFiles = unique(entryFiles);

% Compute dependency closure for each entry file
allDeps = strings(0,1);
missing = strings(0,1);

% Put project on path so requiredFilesAndProducts can resolve dependencies
addpath(genpath(tmpDir));

for i = 1:numel(entryFiles)
    ef = char(entryFiles(i));
    try
        req = matlab.codetools.requiredFilesAndProducts(ef);
        req = string(req);
        allDeps = [allDeps; req]; %#ok<AGROW>
    catch ME
        missing(end+1) = string(ef) + " :: " + string(ME.message); %#ok<AGROW>
    end
end

allDeps = unique(allDeps);

% Keep only deps that are inside the unzipped project (exclude MATLAB/toolbox)
inProject = startsWith(allDeps, string(tmpDir));
projectDeps = allDeps(inProject);

% Copy into outFolder preserving hierarchy
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

copied = strings(0,1);
for i = 1:numel(projectDeps)
    src = char(projectDeps(i));
    rel = erase(src, [tmpDir filesep]);
    dst = fullfile(outFolder, rel);
    dstDir = fileparts(dst);
    if ~exist(dstDir, 'dir')
        mkdir(dstDir);
    end
    copyfile(src, dst);
    copied(end+1) = string(rel); %#ok<AGROW>
end

% Write manifest
manifest = fullfile(outFolder, sprintf('MANIFEST_mainPipeline_%s.txt', mode));
fid = fopen(manifest, 'w');
fprintf(fid, 'Mode: %s\n', mode);
fprintf(fid, 'Source zip: %s\n', zipPath);
fprintf(fid, 'Files copied: %d\n\n', numel(copied));
for i = 1:numel(copied)
    fprintf(fid, '%s\n', copied(i));
end
if ~isempty(missing)
    fprintf(fid, '\n--- Dependency resolution warnings/errors ---\n');
    for i = 1:numel(missing)
        fprintf(fid, '%s\n', missing(i));
    end
end
fclose(fid);

if verbose
    fprintf('Mode "%s": copied %d project files to %s\n', mode, numel(copied), outFolder);
    fprintf('Manifest: %s\n', manifest);
    if ~isempty(missing)
        fprintf('Warning: some entry points could not be analyzed (see manifest)\n');
    end
end

% Best-effort cleanup: leave tmpDir if you want to inspect; comment out if desired
try
    rmpath(genpath(tmpDir));
    rmdir(tmpDir, 's');
catch
end
end
