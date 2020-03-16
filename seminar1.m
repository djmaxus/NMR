function [] = seminar1()
    %% load dataset
    load 'data/misformax';
    dataset = m;
    %% normalize dataset
    mu = mean(dataset,2);
    dataset = dataset - mu;
    val_max = max(dataset,[],'all');
    val_min = min(dataset,[],'all');
    dataset = (dataset - val_min) / (val_max - val_min);
    %% draw example
    f = figure;
    f.WindowState = 'maximized';
    train_index = randi([1,size(dataset,2)],1,1);
    drawImage(dataset(:,train_index),f,'jet');
    title('A field sample from dataset');
    waitforbuttonpress;
    %% train/test split
    uiwait(msgbox('Split to train/test','','modal'));
    data = splitTrainTest(dataset,0.1);
    close(f);
    f = figure;
    f.WindowState = 'maximized';
    
    subplot(1,2,1);
    train_index = randi([1,size(data.train,2)],1,1);
    drawImage(data.train(:,train_index),f,'jet');
    title(['Train sample: ',num2str(train_index),' of ',num2str(size(data.train,2))]);
    
    subplot(1,2,2);
    test_index = randi([1,size(data.test,2)],1,1);
    drawImage(data.train(:,test_index),f,'jet');
    title(['Test sample: ',num2str(test_index),' of ',num2str(size(data.test,2))]);
    waitforbuttonpress;
    %% POD
    uiwait(msgbox('Calculate POD','','modal'));
    POD = makePOD(data.train,false,true,0.99);
    %% set up random masks
    uiwait(msgbox('Mask test data','','modal'));
    mask.frac = 0.1;
    mask.num = round(mask.frac * size(data.test,1));
    mask.test = makeMask(data.test,mask.num);
    mask.P.test = cell(size(mask.test,2),1);
    for i = 1:size(mask.test,2)
        mask.P.test{i} = maskMatrix(mask.test(:,i));
    end
    mask.label = 'mask';
    data.masked.test = data.test .* mask.test;
    data.masked.label = 'gappy';
    test_index = randi([1,size(data.test,2)],1,1);
    compPlot(data.test,mask,data.masked,test_index);
    waitforbuttonpress;
    %% oblique projection
    uiwait(msgbox('Calculate oblique projections','','modal'));
    POD.O.test = cell(size(data.test,2),1);
    wb = waitbar(0,'Calculating oblique projections');
    for i = 1:size(data.test,2)
        waitbar(i/size(data.test,2),wb);
        POD.O.test{i} = oblique(mask.P.test{i},POD.U(:,1:POD.N));
    end
    waitfor(wb);
    %% restore fields
    wb = waitbar(0,'Reconstruction');
    data.restored.test = zeros(size(data.test));
    for i = 1:size(data.test,2)
        waitbar(i/size(data.test,2),wb);
        data.restored.test(:,i) = POD.O.test{i} * data.test(:,i);
    end
    data.restored.label = 'reconstructed';
    waitbar(1,wb,'Close and look what we got!');
    waitfor(wb);
    
    compPlot(data.test,data.masked,data.restored,test_index);
    waitforbuttonpress;
    %% greedy masks
    uiwait(msgbox('Want a better mask?','','modal'));
    grd.num = POD.N;
    grd.frac = grd.num/size(data.test,1);
    [~,grd.P.test] = greedy(POD.U(:,1:POD.N),POD.N);
    grd.test = sum(grd.P.test,2);
    data.masked_grd.test = data.test .* grd.test;
    data.masked_grd.label = 'gappy + greedy';
    POD.Ogrd.test = oblique(grd.P.test,POD.U(:,1:POD.N));
    
    wb = waitbar(0,'Reconstruction');
    data.res_grd.test = zeros(size(data.test));
    for i = 1:size(data.test,2)
        waitbar(i/size(data.test,2),wb);
        data.res_grd.test(:,i) = POD.Ogrd.test * data.test(:,i);
    end
    data.res_grd.label = 'reconstructed (greedy)';
    waitbar(1,wb,'Close and look what we got!');
    waitfor(wb);
    compPlot(data.test,data.masked_grd,data.res_grd,test_index);
    %% draw MSE
    % plotMSE(data,POD,mask,grd);
    
    
end
