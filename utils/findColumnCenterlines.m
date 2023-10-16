function [skeleton,skeletonWithClosing]=findColumnCenterlines(fftImage)
%% Display the original gray scale image.
fftImage=mat2gray(fftImage);
subplot(1, 3, 1);
imshow(fftImage, []);
title('Raw FFT image');

% Now we're ready to begin, since we have our binary image.
%% Skeletonized and thinned to remove spurious branches
skeleton = bwmorph(fftImage,'thin',inf);

% Display it.
subplot(1, 3, 2);
imshow(skeleton, []);
title('Column centers (skeleton)');

%optional, will seal any gaps but adds pixels that may not have been there
morphClose = bwmorph(skeleton,'close');
skeletonWithClosing = bwmorph(morphClose,'thin');
% Display it.
subplot(1, 3, 3);
imshow(skeletonWithClosing, []);
title('Column centers (skeleton + morph. closing)');

%% For distance from center
%{
% Now get the distance transform.
edtImage = bwdist(~fftImage);
%edtImage=edtImage>prctile(edtImage,0);
% Display it.
subplot(2, 2, 3);
imshow(edtImage, []);
title('Euclidean Distance Transform', 'FontSize', fontSize);
% Now multiply
distanceFromEdge = edtImage .* single(skeletonizedImage);
% Display it.
subplot(2, 2, 4);
imshow(distanceFromEdge, []);
title('Distance From Edge', 'FontSize', fontSize);
%}
end