function [] = drawImage(ax,data,clrmap)
    %DRAWIMAGE Summary of this function goes here
    %   Detailed explanation goes here
    imagesc(ax,reshape(data,60,60));
    colormap(ax,clrmap);
    colorbar(ax);
    axis(ax,'square');
end

