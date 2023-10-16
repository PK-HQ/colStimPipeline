function imgsc(varargin)
img=double(varargin{1});
if size(varargin,2)>1
    titleStr=varargin{2};
else
    titleStr='';
end


%% imgsc(img, cmapStr, title)
nanColor=1;
cmin=min(img(:),[],'omitnan');
cmax=max(img(:),[],'omitnan');


%grey nan
imAlpha=ones(size(img));
imAlpha(isnan(img))=0;

if size(imAlpha,3)>1
    imAlpha=mean(imAlpha,3)>0;
end
imagesc(img,'AlphaData',imAlpha);
%caxis([cmin cmax]);
set(gca,'color',nanColor*[1 1 1]);
%
%imagesc(img,[cmin cmax])
if size(img,1) == size(img,2)
    axis square
else
    axis image
end
%colorbar
[colormapRB,colormapR,colormapB] = fireice;

if cmax<0
    colormap(colormapB)
elseif cmin>=0
    colormap(colormapR)
elseif cmin<0 & cmax>0
    colormap(colormapRB)
else
    colormap(gray)
end

colormap(gray)
if contains(titleStr,'\') || contains(titleStr,'{')
    title(titleStr,'FontWeight','Normal','Interpreter','tex')
else
    title(titleStr,'FontWeight','Normal')
end
set(gcf,'color','w')
end