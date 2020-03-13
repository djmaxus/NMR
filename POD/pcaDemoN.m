function [fig] = pcaDemoN(dataset,varargin)
    %PCADEMON Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
        fig = varargin{1};
    else
        fig = figure;
    end
    fig.Name = 'N-D data';
    xmin = min(dataset,[],'all');
    xmax = max(dataset,[],'all');
    % train/test split
    data = splitTrainTest(dataset,0.1);

    subplot(1,2,1);
    train_index = randi([1,size(data.train,2)],1,1);
    drawImage(data.train(:,train_index),fig,'jet');
    caxis([xmin xmax]);
    title(['Train sample: ',num2str(train_index),' of ',num2str(size(data.train,2))]);

    subplot(1,2,2);
    test_index = randi([1,size(data.test,2)],1,1);
    drawImage(data.train(:,test_index),fig,'jet');
    caxis([xmin xmax]);
    title(['Test sample: ',num2str(test_index),' of ',num2str(size(data.test,2))]);
end

