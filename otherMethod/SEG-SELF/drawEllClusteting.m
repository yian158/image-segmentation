%Draws an Ellipse set 
function [hfig] = drawEllClusteting(Iclust,EL,dy,dx)
maxClust = max(max(Iclust));
cmap = hsv(maxClust+1);
cmap(1,:) = [1 1 1];

hfig = figure;
image(Iclust+1);
colormap(cmap);
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