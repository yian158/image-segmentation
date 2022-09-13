%Reads an image - resize - coverts it to binary 
function [I] = myImRead(fname,per) 

I = imread(fname);
if per ~= 1,
    I = imresize(I, per, 'nearest');
end

c(1) = median(double(I(1,:)));
c(2) = median(double(I(:,1)));
c(3) = median(double(I(:,size(I,2))));

bc = median(c);

for i=1:size(I,1),
    for j=1:size(I,2),
        if I(i,j) ~= bc,
            I(i,j) = 1;
        else
            I(i,j) = 0;
        end
    end
end
