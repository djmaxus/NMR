function [fig] = drawImage(data,fig,clrmap)
    %DRAWIMAGE Summary of this function goes here
    %   Detailed explanation goes here
    fig = imagesc(reshape(data,60,60));
    fig = colormap(gca,clrmap);
    colorbar;
    axis square;
end

