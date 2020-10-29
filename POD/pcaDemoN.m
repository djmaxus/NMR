function [] = pcaDemoN(fig)
    %PCADEMON Summary of this function goes here
    %   Detailed explanation goes here
    clf(fig);
    ax = axes(fig);
    %%
    % train/test split
    [data.train,data.test,~] = genData(60,1);
    xmin = min(data.train,[],'all');
    xmax = max(data.train,[],'all');
    
    subplot(1,2,1);
    train_index = randi([1,size(data.train,2)],1,1);
    drawImage(ax,data.train(:,train_index),'jet');
    caxis([xmin xmax]);
    title(['Train sample: ',num2str(train_index),' of ',num2str(size(data.train,2))]);

    subplot(1,2,2);
    test_index = randi([1,size(data.test,2)],1,1);
    drawImage(data.train(:,test_index),fig,'jet');
    caxis([xmin xmax]);
    title(['Test sample: ',num2str(test_index),' of ',num2str(size(data.test,2))]);
end

