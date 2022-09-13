%testSeg: only check the segmentation stage 
function [ IClustTotal,totEll,INITSEG] = runMainAlgo(Iorig,AICBIC_SELECTION,testSeg,METHODSEG,NeighborhoodSize )
I = Iorig;
I = imresize(I, 0.5);

[INITSEG,L] = getSegmentation(I,METHODSEG,NeighborhoodSize);

if testSeg == 0, %check only the segmentation stage 
    INITSEG = imresize(INITSEG, [size(Iorig,1),size(Iorig,2)], 'nearest');
    IClustTotal = imresize(double(L), [size(Iorig,1),size(Iorig,2)], 'nearest');
    totEll = [];
    return;
end

S = regionprops(L,'Area');
for i=1:length(S),
    area(i) = S(i).Area;
end
MArea = mean(area(area > 250));
IClustTotal = zeros(size(I,1),size(I,2));
tStart = tic;
for i=1:max(L(:)),
    O = (L == i);
    
    s = regionprops(O,'BoundingBox');
    apoX = round(s.BoundingBox(2));
    eosX = min(size(O,1),apoX+round(s.BoundingBox(4)));
    apoY = round(s.BoundingBox(1));
    eosY = min(size(O,2),apoY+round(s.BoundingBox(3)));
    Ocrop = O(apoX:eosX,apoY:eosY);
    
    if area(i) < MArea,
        NUMEllipses = 1;
        IClust = Ocrop;
        totEll(i).EL = [];
    else
        [IClust,EL,NUMEllipses] = runMergeFitting(Ocrop,AICBIC_SELECTION);%DEFA method
        totEll(i).EL = EL;
        totEll(i).NUMEllipses = NUMEllipses;
    end
    M = max(IClustTotal(:));
    Bit = IClust > 0;
    IClust = IClust+M*Bit;
    IClustTotal(apoX:eosX,apoY:eosY) = IClustTotal(apoX:eosX,apoY:eosY)+IClust;
end
INITSEG = imresize(INITSEG, [size(Iorig,1),size(Iorig,2)], 'nearest');
IClustTotal = imresize(IClustTotal,[size(Iorig,1),size(Iorig,2)],'nearest');
tElapsed = toc(tStart)
end

