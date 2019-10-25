function [] = compPlot(data,masked,restored,index)
    %COMPPLOT Summary of this function goes here
    %   Detailed explanation goes here
    sample_name = ['Sample ',num2str(index)];
    fig = figure;
    fig.WindowState = 'maximized';
    
    subplot(1,3,1);
    drawImage(data(:,index),fig);
    title([sample_name,': true image']);
    
    subplot(1,3,2);
    drawImage(masked.test(:,index),fig);
    title([sample_name,': masked']);
    
    subplot(1,3,3);
    drawImage(restored.test(:,index),fig);
    title([sample_name,': reconstructed']);
end

