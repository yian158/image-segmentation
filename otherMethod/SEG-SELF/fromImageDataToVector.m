% convert imageData of specific value to vector 
%V(:,1): x cordinates
%V(:,2): y cordinates

function [Points] = fromImageDataToVector(Im,val)

[x,y] = find(Im == val);
N = length(x);
Points = zeros(N,2);
Points(:,1) = x;
Points(:,2) = y;


