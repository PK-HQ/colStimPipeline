function vasculatureImgSegmented=vascSegmentation(vasculatureImg,method)
%% Segments vasculature image for coregistration, choice of methods

% convert to double if not
if ~isa(vasculatureImg,'double')
    vasculatureImg=double(vasculatureImg);
end
    
switch method
    case {'frangi'}
        %% Vascular segmentation
        % ----- Mask positive vasculature -----
        %       .FrangiScaleRange : The range of sigmas used, default [1 8]
        %       .FrangiScaleRatio : Step size between sigmas, default 2
        %       .FrangiBetaOne : Frangi correction constant, default 0.5
        %       .FrangiBetaTwo : Frangi correction constant, default 15
        %       .BlackWhite : Detect black ridges (default) set to true, for
        %                       white ridges set to false.
        options.FrangiScaleRange=[2 12]; %higher range emphasizes large vessels
        options.FrangiScaleRatio=.5;
        options.FrangiBetaOne=2; %may not do much
        options.FrangiBetaTwo=7; %higher will reduce small vessels
        options.BlackWhite=true;
        [vasculatureImgSegmented] = FrangiFilter2D(vasculatureImg,options);
end
end