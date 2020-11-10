close all;
clear;
%% load dataset
load 'data/misformax2';
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
drawImage(dataset(:,1),f);
title('A field sample from dataset');
waitforbuttonpress;
%% train/test split
uiwait(msgbox('Split to train/test','','modal'));
data = splitTrainTest(dataset,0.1);
%% POD
uiwait(msgbox('Calculate POD','','modal'));
POD = makePOD(data.train,0.9);
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
compPlot(data.test,mask,data.masked,1);
waitforbuttonpress;
%% oblique projection
uiwait(msgbox('Calculate oblique projections','','modal'));
POD.O.test = cell(size(data.test,2),1);
wb = waitbar(0,'Calculating oblique projections');
for i = 1:size(data.test,2)
    waitbar(i/size(data.test,2),wb);
    POD.O.test{i} = oblique(mask.P.test{i},POD.Ut);
end
%% restore fields
uiwait(msgbox('Reconstruct data','','modal'));
waitbar(0,wb,'Reconstruction');
data.restored.test = zeros(size(data.test));
for i = 1:size(data.test,2)
    waitbar(i/size(data.test,2),wb);
    data.restored.test(:,i) = POD.O.test{i} * data.test(:,i);
end
data.restored.label = 'reconstructed';
waitbar(1,wb,'Close and look what we got!');
waitfor(wb);
compPlot(data.test,data.masked,data.restored,1);
waitforbuttonpress;
%% greedy masks
uiwait(msgbox('Want a better mask?','','modal'));
grd.num = size(POD.S,2);
grd.frac = grd.num/size(data.test,1);
grd.P.test = greedy(POD.U(:,1:grd.num));
grd.test = sum(grd.P.test,2);
data.masked_grd.test = data.test .* grd.test;
data.masked_grd.label = 'gappy + greedy';
POD.Ogrd.test = oblique(grd.P.test,POD.Ut);

wb = waitbar(0,'Reconstruction');
data.res_grd.test = zeros(size(data.test));
for i = 1:size(data.test,2)
    waitbar(i/size(data.test,2),wb);
    data.res_grd.test(:,i) = POD.Ogrd.test * data.test(:,i);
end
data.res_grd.label = 'reconstructed (greedy)';
waitbar(1,wb,'Close and look what we got!');
waitfor(wb);
compPlot(data.test,data.masked_grd,data.res_grd,1);
%% draw MSE
% plotMSE(data,POD,mask,grd);


