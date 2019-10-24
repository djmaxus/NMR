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
pause(1);
%% train/test split
data = splitTrainTest(dataset,0.1);
%% POD
POD = makePOD(data.train,0.8);
%% draw MSE
plotMSE(data,POD);