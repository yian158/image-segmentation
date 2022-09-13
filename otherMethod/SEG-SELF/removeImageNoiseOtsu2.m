function [ A,L ] = removeImageNoiseOtsu2( I,toPlot )

[ A0,L0,CC,level ] = removeImageNoise( I,toPlot );
S = regionprops(CC, 'Area');
for i=1:length(S),
    v(i) = S(i).Area;
end

T1 = sqrt(var(v)) / mean(v)
M = max(I(:));
I2 = I;
I2(I > level*M) = level*M;

[ A1,L1,CC1,level1 ] = removeImageNoise( I2,toPlot );
v = [];
S = regionprops(CC1, 'Area');
for i=1:length(S),
    v(i) = S(i).Area;
end

T2 = sqrt(var(v)) / mean(v)

A = A0;
L = L0;

if T2 < T1,
    A = A1;
    L = L1;
end

if toPlot == 1,
    figure;
    imagesc(A0);
    title('Otsu 1');
    figure;
    imagesc(A1);
    title('Otsu 2');
    
end

