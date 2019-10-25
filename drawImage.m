function [] = drawImage(data,fig)
    %DRAWIMAGE Summary of this function goes here
    %   Detailed explanation goes here
    fig = imagesc(reshape(data,60,60));
    colormap 'jet';
    colorbar;
    caxis([0 1]);
    axis square;
end

