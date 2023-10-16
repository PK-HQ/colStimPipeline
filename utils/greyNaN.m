function greyNaN(varargin)
img=varargin{1};
if size(varargin,2)>1
    nanColor=varargin{2};
else
    nanColor=1;
end

if isequal(class(img),'logical')
    img=double(img);img(img==0)=NaN;
elseif isequal(class(img),'uint8')
    img=double(img);%img(img==0)=NaN;
elseif sum(isnan(img))>0
    %nothing
end

imAlpha=ones(size(img));
imAlpha(isnan(img))=0;

if size(imAlpha,3)>1
    imAlpha=mean(imAlpha,3)>0;
end
imagesc(img,'AlphaData',imAlpha);
set(gca,'color',nanColor*[1 1 1]);
end