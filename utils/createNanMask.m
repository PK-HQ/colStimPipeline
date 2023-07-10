function nanMask=createNanMask(Mask)
nanMask=double(Mask); % convert binary mask to NaN mask
nanMask(nanMask==0)=NaN;
end