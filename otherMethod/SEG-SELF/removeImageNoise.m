function [ A,L ] = removeImageNoise( I,toPlot )
toPlot = 0;
level = graythresh(I); %Otsu Method
disp('Level Otsu:');
disp(level);
A = imbinarize(I,level);

CC = bwconncomp(A, 8);
S = regionprops(CC, 'Area');
for i=1:length(S),
    v(i) = S(i).Area;
end
v = sort(v,'descend');
v = v(v > 20);
TreshP = 0.1*median(v);
L = labelmatrix(CC);
A = ismember(L, find([S.Area] >= TreshP));
A = imfill(A,'holes');

CC = bwconncomp(A, 8);
L = labelmatrix(CC);

if toPlot ==1,
    figure;
    imagesc(I);
    title('Original Image');
    
    figure;
    imagesc(A);
    title('After Denosing');
end