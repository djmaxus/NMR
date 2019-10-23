%% load dataset
% load 'data/gauss60_5+';
load 'data/misformax';
dataset = m;
%% normalize dataset
mu = mean(dataset,2);
dataset = dataset - mu;

%% split dataset
i_test = randi([1,size(dataset,2)],10,1);
i_test = sort(unique(i_test));
test = zeros(size(dataset,1),length(i_test));
train = dataset;
for i = 1:length(i_test)
    test(:,i) = dataset(:,i_test(i));
    train(:,i_test(i)-i) = [];
end
%% perform SVD
[U,S,V] = svd(train);
%% mean squared error
figure;
subplot(1,2,1);
semilogy(MSE(U,S,train));
hold on;
semilogy(MSE(U,S,test));

%% spectrum
diagS = diag(S);
Nt = eigEnNum(diagS,0.9);
subplot(1,2,2);
semilogy(diagS(1:end-1));
%% set gap parameters
% unmarred fraction
p = 0.1;
% unmarred number
n_unmarred = round(p * size(train,1));

%% set up masks
mask = makeMask(train,n_unmarred);
mask_test = makeMask(test,n_unmarred);
for j = 1:size(train,2)
    M{j} = maskMatrix(mask(:,j));
end
%% marred average
data_masked = train .* mask;
mu_masked = sum(data_masked,2)./ (sum(mask,2) + eq(sum(mask,2),0));
data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* mu_masked;

%% iterations
max_iter = 20;
for num_iter = 1:max_iter
    [Um,Sm,~] = svd(data_masked);
    hold on;
    diagS = diag(Sm);
    f = semilogy(diagS(1:end-1));
    %% least squares problem
    A =  Um(:,1:Nt);
    a = zeros(Nt,size(data_masked,2));
    for j = 1:size(data_masked,2)
        B = data_masked(:,j);
        x = lscov(M{j}*A,M{j}*B);
        a(:,j) = x;
    end
    
    data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* (A*a);
end
figure;
subplot(2,1,1);
out = dataset(:,1);
imagesc(reshape(out,60,60));
subplot(2,1,2);
out = data_masked(:,1);
imagesc(reshape(out,60,60));
%% repair test fields by unmarried train data
A = U(:,1:Nt);
a = zeros(Nt,1);
B = test(:,1);
MM = maskMatrix(mask_test(:,1));

x = lscov(MM*A,MM*B);
a(:,j) = x;
restored = A * x;
figure;
subplot(2,1,1);
out = test(:,1);
imagesc(reshape(out,60,60));
subplot(2,1,2);
out = restored(:,1);
imagesc(reshape(out,60,60));
%% oblique projection
P = MM';
A = U(:,1:Nt);
PP = P'*A;
[Up,Sp,Vp] = svd(PP);
Spp = zeros(size(Sp'));
dSp = diag(Sp);
dSp = eye(length(dSp))./ dSp;
Spp(1:end,1:length(dSp))= dSp;
PPp = Vp * Spp * Up';
ObProj = A * PPp * P';
%% check if entry is non-unique in array
function [flag] = is_repeat(entry_index,array)
sum_left = sum(eq(array(entry_index),array(1:entry_index-1)));
sum_right = sum(eq(array(entry_index),array(entry_index+1:end)));
flag = (sum_left + sum_right) > 0;
end

%% calculate mask
function mask = makeMask(dataset,n_unmarred)
mask = zeros(size(dataset,1),size(dataset,2));
i_unmarred = randi([1,size(dataset,1)],n_unmarred,size(dataset,2));
for i = 1:size(dataset,2)
    for j = 1:n_unmarred
        while is_repeat(j,i_unmarred(:,i))
            val = randi([1,size(dataset,1)],1,1);
            i_unmarred(j,i) = val;
        end
        mask(i_unmarred(j,i),i) = 1;
    end
end
end
function out = maskMatrix(mask)
out = zeros(sum(mask,1),size(mask,1));
pos = 1;
for i = 1:size(mask,1)
    if mask(i) == 1
        out(pos,i) = 1;
        pos = pos + 1;
    end
end
end
function out = eigEnNum(S,fraction)
total = sum(S);
for i = 1:length(S)
    if sum(S(1:i)) >= fraction * total
        out = i;
        break;
    end
end
end
function [out] = MSE(U,S,data)
out = zeros(1,size(S,2));
for n = 1:size(S,2)-1
    A = U(:,1:n);
    restore = A * A' * data;
    delta = data - restore;
    delta = sum(delta.^2,'all')/numel(data);
    out(n) = delta;
end
end
function [M,F] = calcMF(data,U)
    M = U' * U;
    F = U' * data;
end