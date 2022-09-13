function [ok] = drawDistEllClustetingMerge(Dtemp,EL,dy,dx)

ok = figure;
imagesc(Dtemp);
for k=1:numel(EL),
    if EL(k).Label == k,
        theta = -EL(k).phi*pi/180;
        x0 = EL(k).C(1)-dx;
        y0 = EL(k).C(2)-dy;
        a = EL(k).a;
        b = EL(k).b;
        PE = ellipseToPolygon([x0 y0 a b theta]);
        hold on;
        plot(PE(:,1),PE(:,2),'k','LineWidth',2);
        hold on;
        text(x0,y0,sprintf('%d',EL(k).Label));
    end
end