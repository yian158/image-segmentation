%It runs the AEFA (Splitting) method

%I: original binary image
%AICBIC_SELECTION: Set AICBIC_SELECTION = 1, to use AIC is selected else BIC is used
%IClust:Clustering of image pixels
%EL: Ellipses 
%NUMEllipses:Number of Ellipses 

function [IClust,EL,NUMEllipses] = runSlitFitting(I,AICBIC_SELECTION)
toPlot = 0;
Iorig = I;
lines = size(I,1);
cols = size(I,2);
I = zeros(3*lines,3*cols);
I(lines+1:2*lines,cols+1:2*cols) = Iorig;
IClust = I;

[nCompl] = getObjectComplexity(I);

[EL,NUMEllipses,area,AIC,BIC,minAICBIC,SI,hfig1] = runSkeletonInitFAST_DIST(IClust,0,lines,cols,nCompl,AICBIC_SELECTION);

[EL,IClust,~,TotalPerf] = runEllClustering(EL,IClust,area);
if toPlot == 1,
    [ok] = drawEllClusteting(IClust(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
    title(sprintf('%4.2f',100*TotalPerf));
 
    
    drawEllClusteting(IClust(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
    for i=1:NUMEllipses,
        hold on;
        text(EL(i).C(1) -cols,EL(i).C(2)-lines,sprintf('%d',i));
    end
    title(sprintf('%4.2f',100*TotalPerf));
    
    figure;
    plot(AIC,'-.o');
    hold on;
    plot(BIC,'--rs');
    legend('AIC','BIC');
end
perf = [EL.outPixels];
[min_perf, pos] = max(perf);
IClust = IClust(lines+1:2*lines,cols+1:2*cols);







