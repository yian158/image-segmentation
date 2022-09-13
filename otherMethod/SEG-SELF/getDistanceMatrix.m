%Computes the distance matrix of points under the user metric 
%metric: 2 Euclidean Distance
%metric: 3 Inner distance
%metric: 4 Modified Inner distance
%metric: 5 Line distance

function [Dist] = getDistanceMatrix(Points, Im,val,metric)

if metric == 2,
    [Dist] = getEuclideanDistanceMatrix(Points);
elseif metric == 3,
    [Dist] = getInnerDistanceMatrix2DData(Points);
elseif metric == 4,
    [Dist_E] = getEuclideanDistanceMatrix(Points);
    [Dist_I] = getInnerDistanceMatrix2DData(Points);
    N = size(Points,1);
    Dist = zeros(N,N);
    for i=1:N,
        for j=i+1:N,
            r = Dist_I(i,j) / Dist_E(i,j);
            Dist(i,j) = Dist_I(i,j)^r;
            Dist(j,i) = Dist(i,j);
        end
    end
elseif metric == 5,
    [Dist] = getLineInnerDistanceMatrix2DData(Points,Im,val);
end


function [Dist] = getEuclideanDistanceMatrix(Points)

N = size(Points,1);
Dist = zeros(N,N);

for i=1:N,
    p1 = Points(i,:);
    for j=i+1:N,
        p2 = Points(j,:);
        Dist(i,j) = norm(p1-p2);
        Dist(j,i) = Dist(i,j);
    end
end


% Computes the  InnerDistanceMatrix under 2D data (8-connections)
%

function [Dist] = getInnerDistanceMatrix2DData(Points)


N = size(Points,1);
Dist1 = zeros(N,N);
X = Points(:,1);
Y = Points(:,2);

for i=1:N,
    x = Points(i,1);
    y = Points(i,2);
    for i1=-1:1,
        for i2=-1:1,
            if i1 ~= 0 && i2 ~= 0
                x1 = x+i1;
                y1 = y+i2;
                pos = find(X == x1 & Y == y1);
                if ~isempty(pos),
                    Dist1(i,pos) = norm([i1 i2]);
                    Dist1(pos,i) = Dist1(i,pos);
                end
            end
        end
    end
end

[Dist] = graphallshortestpaths(sparse(Dist1));


function [Dist] = getLineInnerDistanceMatrix2DData(Points,Im,val)

N = size(Points,1);
Dist = zeros(N,N);
N = size(Points,1);

for i=1:N,
    p1 = Points(i,:);
    
    for j=i+1:N,
        p2 = Points(j,:);
        dp = p2-p1;
        dist = norm(dp);
        numPoints = ceil(dist+1);
        
        for k=1:numPoints-1,
            p = round(p1 + (k/numPoints)*dp);
            if Im(p(1),p(2)) ~= val,
                dist = 10^20; 
                break;
            end
        end
        Dist(i,j) = dist;
        Dist(j,i) = dist;
    end
end


