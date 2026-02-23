function score = registrationQuality(fixed, movingReg)
    fixed = im2double(fixed);  movingReg = im2double(movingReg);

    % Intensity
    score = ssim(movingReg,fixed);
    % gradient correlation ( blood vessels)
   % score.gradcorr = corr2(imgradient(fixed), imgradient(movingReg));

end
