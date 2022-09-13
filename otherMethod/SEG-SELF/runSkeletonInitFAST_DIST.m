%Initiliazation of Ellipses based on Skeleton of A and distance tranform of
%A
function [EL,NUMEllipses,area,AIC,BIC,minAICBIC,SI,hfig1] = runSkeletonInitFAST_DIST(A,toPlot,lines,cols,nCompl,AICBIC_SELECTION)
Cent = [];
Rad = [];
hfig1 = [];
%SK = bwmorph(A,'skel',inf);
SK = bwmorph(A,'thin',inf);%????
EP = bwmorph(SK,'endpoints');

BD = SK.*bwdist(1-A);
FULLPLOT = 0;
stats = regionprops(double(A), 'Area');
area = stats.Area;

[Points(:,1) Points(:,2)] = find(BD > 0);
W = zeros(1,size(Points,1));
for i=1:length(W),
    W(i) = BD(Points(i,1),Points(i,2));
end

%LOCAL MAXIMA COMPUTATION 
LM = zeros(1,size(Points,1));
VAL = zeros(1,size(Points,1));
ite = 0;
for i=1:size(Points,1),
    v = [(Points(:,1)-Points(i,1)).^2+(Points(:,2)-Points(i,2)).^2];
    maxDist = sqrt(max(v)/(mean(v)+0.001));
    
    ite = ite+1;
    LM(ite) = i;
    VAL(ite) = BD(Points(i,1),Points(i,2))+(1/(1+1000*maxDist));
end
LM = LM(1:ite);
VAL = VAL(1:ite);

[SVAL,pos] = sort(VAL,'descend');
ite = 1;
i = 1;
Cent(ite,1:2) = Points(LM(pos(i)),1:2);
Rad(ite) = SVAL(i);
%Selection of candicate circle centers 
for i=2:length(VAL),
    v = sqrt([(Cent(:,1)-Points(LM(pos(i)),1)).^2+(Cent(:,2)-Points(LM(pos(i)),2)).^2]);
    if SVAL(i) < 0.03*SVAL(1)+1,
        break;
    end
    ok = 1;
    x0 = Points(LM(pos(i)),1);
    y0 = Points(LM(pos(i)),2);

    for j=1:length(v),
         mdist = v(j);
        if EP(x0,y0) == 0 && mdist <= Rad(j)+SVAL(i),%min(Rad(j),SVAL(i))+1%Rad(j)+SVAL(i)+1, %out of cirles; %max(Rad(j),SVAL(i))+1,  %
            ok = 0;
            break;
        elseif EP(x0,y0) == 1 && mdist <= Rad(j),
            ok = 0;
            break;
        end
    end
    if ok == 1,
        ite = ite+1;
        Cent(ite,1:2) = Points(LM(pos(i)),1:2);
        Rad(ite) = SVAL(i);
    end

end
%Selection of  circle centers based on AIC
[EL,~] = initEll(Cent,Rad,1);
[EL,~,~,~] = runEllClustering(EL,A,area);
[EL,IClust,BestDTemp,TotalPerf] = runEllClustering(EL,A,area);
%TotalPerf = sum([EL.InArea]) / area;

[AIC(1),BIC(1),RES_AICBIC,bestAICBIC,SI(1,1:3)] = getAIC_BIC(nCompl,TotalPerf,1,AICBIC_SELECTION,IClust,EL);

if FULLPLOT == 1,
    [~] = drawDistEllClusteting(BestDTemp(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
end
minAICBIC = RES_AICBIC; 
ELLSET = 1;
BEST_ELLSET = ELLSET;
ite = 1;
CANDICATE_SET = [2:length(Rad)];
BETTERSOL = 0;

while ~isempty(setdiff(CANDICATE_SET,ELLSET)),

    CANDICATE_SET = getSortedCANDICATE_SET_EM2(CANDICATE_SET,ELLSET,BestDTemp,Cent,Rad);
    if FULLPLOT == 1,
        [ok] = drawCENT(CANDICATE_SET,ELLSET,EL,Cent,BestDTemp(lines+1:2*lines,cols+1:2*cols),lines,cols);
    end
    change = 0;
    if ~isempty(find(ELLSET == CANDICATE_SET(1), 1)),
        break;
    end
    for id=1:length(CANDICATE_SET),
        val = CANDICATE_SET(id);
       
        ELLSET = [ELLSET val];
        [EL,NUMEllipses] = initEll(Cent,Rad,ELLSET);
        [EL,IClust,DTemp,TotalPerf] = runEllClustering(EL,A,area);      
        %TotalPerf = sum([EL.InArea]) / area;
        [AIC1,BIC1,RES_AICBIC,bestAICBIC,SI1(1,1:3)] = getAIC_BIC(nCompl,TotalPerf,NUMEllipses,AICBIC_SELECTION,IClust,EL);
        
        %Centers correction (CC) 
        [newELLSET] = getNEWELLSET(Cent,Rad,ELLSET,CANDICATE_SET,IClust,DTemp);
        if sum(abs(ELLSET-newELLSET)) > 0,
            disp('Changing in ELLSET');
            ELLSET
            newELLSET
            [newEL,~] = initEll(Cent,Rad,newELLSET);
            [newEL,IClustTemp,newDTemp,newTotalPerf] = runEllClustering(newEL,A,area);
            %newTotalPerf = sum([newEL.InArea]) / area;
            newTotalPerf
            TotalPerf
            if newTotalPerf > TotalPerf,
                disp('Better  TotalPerf. Changing done.')
                [AIC1,BIC1,RES_AICBIC,bestAICBIC,SI1(1,1:3)] = getAIC_BIC(nCompl,newTotalPerf,NUMEllipses,AICBIC_SELECTION,IClustTemp,newEL);
                EL = newEL;
                DTemp = newDTemp;
                ELLSET = newELLSET;
                TotalPerf = newTotalPerf;
            end
        end
        
        ite = ite+1;
        if RES_AICBIC < minAICBIC,
            BEST_ELLSET = ELLSET;
            minAICBIC = RES_AICBIC;
            change = 1;
            AIC(length(ELLSET)) = AIC1;
            BIC(length(ELLSET)) = BIC1;
            SI(length(ELLSET),1:3) = SI1(1,1:3);
            BestDTemp = DTemp;
            BestTotalPerf = TotalPerf;
            BestEL = EL;
            break;
        else
            AIC(length(ELLSET)) = AIC1;
            BIC(length(ELLSET)) = BIC1;
            SI(length(ELLSET),1:3) = SI1(1,1:3);
            BETTERSOL = BETTERSOL+1;
            BestDTemp = DTemp;
            %ELLSET = ELLSET(1:NUMEllipses-1);
            break;
        end 
    end
   %bestAICBIC
   %minAICBIC
  
    if change == 1,
        BETTERSOL = 0;
        ELLSET = BEST_ELLSET;
        CANDICATE_SET = setdiff(CANDICATE_SET,BEST_ELLSET);
         if FULLPLOT == 1,
            [~] = drawDistEllClusteting(BestDTemp(lines+1:2*lines,cols+1:2*cols),BestEL,lines,cols);
            title(sprintf('TotalPerf = %2.4f - %d',BestTotalPerf,length(BEST_ELLSET)));
         end
    elseif BETTERSOL <= 12 && ~isempty(CANDICATE_SET) && bestAICBIC < minAICBIC,%?????10
        if FULLPLOT == 1,
            [~] = drawDistEllClusteting(BestDTemp(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
            title(sprintf('(Retry %d), TotalPerf = %2.4f - %d',BETTERSOL,TotalPerf,length(TotalPerf)));
        end
    else
        BETTERSOL
        CANDICATE_SET
        bestAICBIC
        minAICBIC
        break;
    end
%     if length(ELLSET) == 6,
%             break;
%      end
end

ELLSET = BEST_ELLSET;

[EL,NUMEllipses] = initEll(Cent,Rad,ELLSET);


if toPlot == 1,
    cmap = colormap(jet(256));
    cmap(1,:) = [1 1 1];
    hfig1 = figure;
    imagesc(SK(lines+1:2*lines,cols+1:2*cols));
    hold on;
    imagesc(BD(lines+1:2*lines,cols+1:2*cols)+double(A(lines+1:2*lines,cols+1:2*cols)));
    colormap(cmap);
    title(sprintf('Distances - Complexity = %2.2f',nCompl));
    
    for i=1:length(ELLSET),
            hold on;
            text(Cent(ELLSET(i),2)-cols,Cent(ELLSET(i),1)-lines,sprintf('%d',i));
    end
    cid = length(ELLSET)+1;
    for i=1:size(Cent,1),
        if isempty(find(ELLSET == i, 1)),
            hold on;
            text(Cent(i,2)-cols,Cent(i,1)-lines,sprintf('%d',cid));
            cid = cid+1;
        end
    end
    [ok] = drawEllClustetingNoFig(EL,lines,cols);
end



%Initilization of Ellipses 
function [EL,NUMEllipses] = initEll(Cent,Rad,ELLSET)

NUMEllipses = length(ELLSET);
EL = [];
for val=1:length(ELLSET),
    id = ELLSET(val);
    EL(val).a = Rad(id);
    EL(val).b = Rad(id);
    if val == 1,
        EL(val).a = Rad(id);
        EL(val).b = Rad(id);
    end
    EL(val).C(1) = Cent(id,2);
    EL(val).C(2) = Cent(id,1);
    EL(val).phi = 0;
    EL(val).InArea = 0;
    EL(val).outPixels = 0;
    EL(val).tomh_area = 0;
    EL(val).tomh_enwsh = 0;
    EL(val).Label = val;
    EL(val).ELLSET = id;
end

function [ok] = drawCENT(CANDICATE_SET,ELLSET,EL,Cent,Dtemp,dy,dx)
figure;
imagesc(Dtemp);
[ok] = drawEllClustetingNoFig(EL,dy,dx);
for id=1:size(Cent,1),
    hold on;
    i = CANDICATE_SET(id);
    if ~isempty(find(ELLSET == i)),
        text(Cent(i,2)-dx,Cent(i,1)-dy,sprintf('%d',id),'FontSize',18,'Color',[1 1 1]);
    else
        text(Cent(i,2)-dx,Cent(i,1)-dy,sprintf('%d',id));
    end
end



function [ok] = drawEllClustetingNoFig(EL,dy,dx)
ok = 1;


for k=1:numel(EL),
    hold on;
    x0 = EL(k).C(1)-dx;
    y0 = EL(k).C(2)-dy;
    a = EL(k).a;
    b = EL(k).b;
    theta = -EL(k).phi*pi/180;
    
    PE = ellipseToPolygon([x0 y0 a b theta]);
    plot(PE(:,1),PE(:,2),'k');
end





function [CS] = getSortedCANDICATE_SET_EM2(CANDICATE_SET,ELLSET,BestDTemp,Cent,Rad)
%find the non-covering points
[posX,posY] = find(BestDTemp > 1);
VAL = zeros(1,length(Rad));
d = zeros(1,length(ELLSET));


for id=1:length(CANDICATE_SET),
    cid = CANDICATE_SET(id);
    for i=1:length(posX),
        x = posX(i);
        y = posY(i);
        for j=1:length(ELLSET),
            OAdist = norm([x y] - Cent(ELLSET(j),1:2));
            OXdist = Rad(ELLSET(j));
            d(j) = OAdist/max(OXdist,0.00001);
        end
        OAdist = norm([x y] - Cent(cid,1:2));
        OXdist = Rad(cid);
        dC = OAdist/max(OXdist,0.000001);
        dC_ELLSET = min(d);
        if dC < dC_ELLSET,
            pr1 = 1-((1/BestDTemp(x,y)^2));
            %pr2 = min(1,1/(max(dC,0.0000001)));
            %vec = [1/BestDTemp(x,y) 1/dC];
            
            pr2 = min(1,1/(dC^2));
            
            VAL(cid) = VAL(cid)+pr1*pr2;
        end
    end
end
[~,CS] = sort(VAL,'descend');


%Computes the newELLSET using bitarate matching between the centers and the
%current clustering 
function  [newELLSET] = getNEWELLSET(Cent,Rad,ELLSET,CANDICATE_SET,IClust,BD)
N = length(CANDICATE_SET);
E = length(ELLSET);
C = zeros(E,N);
for i=1:E,
    [px,py] = find(IClust == i);    
    for j=1:N,
        jid = CANDICATE_SET(j);
        v = sqrt([(Cent(jid,1)-px).^2+(Cent(jid,2)-py).^2]);
        pos = find(v <= Rad(jid));
        for k=1:length(pos),
            s = 1-(BD(px(pos(k)),py(pos(k))))^2;
            if s > 0,
                C(i,j) = C(i,j)+s;
            end
        end
    end
end

[c,cost] = munkres(-C);
newELLSET = CANDICATE_SET(c);





