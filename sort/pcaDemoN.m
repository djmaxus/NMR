function [] = pcaDemoN(fig1,fig2)
    %PCADEMON Summary of this function goes here
    %   Detailed explanation goes here
    clf(fig1);
    ax = axes(fig1);
    sp1 = subplot(1,3,1,ax);
    sp2 = subplot(1,3,2,axes(fig1));
    sp3 = subplot(1,3,3,axes(fig1));
    %%
    % train/test split
    [data.train,data.test,~] = genData(60,1);
    xmin = min(data.train,[],'all');
    xmax = max(data.train,[],'all');

    train_index = randi([1,size(data.train,2)],1,1);
    drawImage(sp1,data.train(:,train_index),'jet');
    caxis(sp1,[xmin xmax]);
    title(sp1,['Train sample: ',num2str(train_index),' of ',num2str(size(data.train,2))]);

    drawImage(sp2,data.test,'jet');
    caxis(sp2,[xmin xmax]);
    title(sp2,'Test sample');

    %% pca demo
    pca = PCA(data.train,0.99);
    pca.plotSummary(fig2,'layout','separate','data',data.train);
    drawImage(sp3,pca.truncate(data.test),'jet');
    caxis(sp3,[xmin xmax]);
    title(sp3,'Test sample: projected');
end
