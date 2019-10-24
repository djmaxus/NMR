function [] = compPlot(F,R,P,index)
    %COMPPLOT Summary of this function goes here
    %   Detailed explanation goes here
    orig = reshape(F(:,index),60,60);
    reco = reshape(R(:,index),60,60);
    mask = reshape(sum(P,2),60,60);
    figure;
    subplot(1,3,1);
    imagesc(orig);
    subplot(1,3,2);
    imagesc(mask .*orig);
    subplot(1,3,3);
    imagesc(reco);
end

