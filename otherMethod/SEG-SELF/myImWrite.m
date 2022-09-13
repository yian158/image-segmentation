function [ cmap ] = myImWrite(I,IClustTotal,GT,ResultsDir,fname,toPlot )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
C = max(IClustTotal(:))+1;
cmap = colorcube(C);
cmap(1,:) = [1 1 1];
S = sum(cmap');
pos = find(S == 0);
if length(pos) == 0,
    cmap = colorcube(C);
    S = sum(cmap');
    pos = find(S == 3);
    cmap(1,:) = [1 1 1];
    cmap(pos,:) = [0 0 0];
end
%imwrite(IClustTotal,cmap,sprintf('%s%s.jpg',ResultsDir,fname));
temp = IClustTotal;
temp(GT == 1) = pos-1;
imwrite(temp,cmap,sprintf('%s%sGT.png',ResultsDir,fname));

if toPlot == 1,
    figure; imagesc(I); colormap('gray');
    title('Original');
     figure; imagesc(IClustTotal); colormap(cmap);
     title('Segmentation');
    figure; imagesc(temp); colormap(cmap);
    title('Segmentation+GT');
end



