function stabcfg = autoaskAboutStabilization(fnTS,TS,iFramesNmlTrial)

%%  stabcfg = askAboutStabilization(fnTS,TS,DAVersion)

%% Ask about stabilization
%{
menuselect = menu_mod('Would you like a run image stabilization?', ...
        'No', ...
        'Yes + Saved Cfg (or Default)', ...
        'Yes + Default Cfg', ...
        'Yes + Custom Cfg' ...
       );
   %}
menuselect=2;
switch (menuselect)
  case 2
    dostab    = 1;
    dodefault = 1;
    dosaved   = 1;
    
  case 3
    dostab    = 1;
    dodefault = 1;
    dosaved   = 0;
    
  case 4
    dostab    = 1;
    dodefault = 0;
    dosaved   = 0;
    
  otherwise
    dostab    = 0;
    dodefault = 0;
    dosaved   = 0;
        
end

% %% Ask about stabilization
% dostab = questdlg('Would you like a run image stabilization?', ...
% 	'Image Stabilization', ...
% 	'Yes','No','Yes');
% dostab = strcmp(dostab,'Yes');

stabcfg   = [];

%% Extra configuration dialogue
if dostab

  [PathName,fnRoot] = fileparts(fnTS);
  fnRoot = regexprep(fnRoot,'TS','');

  % load existing file if exists
  fnStabCfg = fullfile(PathName,sprintf('%sStabCfg.mat',fnRoot));
  if (exist(fnStabCfg,'file') && dosaved)
    disp(['Using Saved StabCfg: ',fnStabCfg,' ...']);
%     loadstabcfg = questdlg('Image stabilization configuration file found. Re-use?', ...
%     'Image Stabilization', ...
%     'Yes','No','Yes');  
%     if (strcmp(loadstabcfg,'Yes'))
      stabcfg = load(fnStabCfg);
      stabcfg.src = fnStabCfg;
      stabcfg.ipu = WideField.ipu.ImageStabilizer(stabcfg.ipu);
%     end
  end
  
  % user config
  if (isempty(stabcfg) || ~isfield(stabcfg,'cfg'))        
    
    if (isfield(TS.Header.Imaging,'CameraName'))
      divlinesX = '[]';
    else
      divlinesX = '[128 129 256 257 384 385]';
    end
    
    stabset = struct( ...
      'spatialHP',     1.5, ...
      'iFrameRef',     iFramesNmlTrial(ceil(numel(iFramesNmlTrial)/2)), ...
      'regwin',        [138 138 375 375], ...
      'divlinesX',     divlinesX, ...
      'iFramesStab',   [], ...
      'SizePxl',       TS.Header.Imaging.SizePxl, ...
      'FrameRate',     TS.Header.Imaging.FrameRate, ...
      'detrend',       '' ...
    );
    
    
    if (~dodefault)
      Answer = ...
        inputdlg({'Motion estimate spatial high pass pre-filter (cyc/mm)', ...
                  'Static reference frame index (0 = single common frame)', ...
                  'Static reference ROI [startx starty endx endy]', ...
                  'CCD division lines', ...
                  'Frames used to estimate stabilization parameters ([] for all)', ...
                  'De-trend ('''' = none, ''quasub'' = quadrat fit)'}, ...
                 'Parameters',1, ...
                 {num2str(stabset.spatialHP), num2str(stabset.iFrameRef), ...
                  mat2str(stabset.regwin), divlinesX, ...
                  mat2str(stabset.iFramesStab), stabset.detrend});
%       drawnow;pause(0.1);

      if isempty(Answer)
        Message = 'Canceled!';
        disp(Message);
        return;
      end

      stabset.spatialHP = str2double(Answer{1});
      stabset.iFrameRef = str2double(Answer{2});
      stabset.regwin      = str2num(Answer{3});
      stabset.divlinesX   = str2num(Answer{4});
      stabset.iFramesStab = str2num(Answer{5});
      stabset.detrend     = Answer{6};

      if isempty(stabset.spatialHP)||~isscalar(stabset.spatialHP)||(stabset.spatialHP<0)|| ...
         isempty(stabset.iFrameRef)||~isscalar(stabset.iFrameRef)||(stabset.iFrameRef<0)||(stabset.iFrameRef>TS.Header.Imaging.nFrame)|| ...
         isempty(stabset.regwin)||stabset.regwin(1)<1||stabset.regwin(2)<1||stabset.regwin(3)>TS.Header.Imaging.FrameWidth||stabset.regwin(4)>TS.Header.Imaging.FrameHeight|| ...
         any(stabset.divlinesX<1) || any(stabset.divlinesX>TS.Header.Imaging.FrameWidth)|| ...
         any(stabset.iFramesStab<1) || any(stabset.iFramesStab>TS.Header.Imaging.nFrame) || ...
         ~ismember(stabset.detrend,{'','quadsub'})
        beep;
        return;
      end
    end

    if isempty(stabset.detrend)
      stabset.detrend = 'meanTC';
    end
    
    stabcfg.fnStabCfg  = fnStabCfg;
    stabcfg.iiFrameRef = stabset.iFrameRef;
    stabcfg.cfg = cell(TS.Header.Index.nBLKTrial,1);
    stabcfg.ipu = WideField.ipu.ImageStabilizer( ...
      'spatialHP',     stabset.spatialHP, ...
      'refframe',      stabset.iFrameRef, ...
      'regwinX',       stabset.regwin([1 3]), ...
      'regwinY',       stabset.regwin([2 4]), ...
      'divlinesX',     stabset.divlinesX, ...
      'regressFrames', stabset.iFramesStab, ...
      'SizePxl',       TS.Header.Imaging.SizePxl, ...
      'FrameRate',     TS.Header.Imaging.FrameRate, ...
      'DBleacher',     stabset.detrend ...
      );
    stabcfg.ipu
  end
end  
