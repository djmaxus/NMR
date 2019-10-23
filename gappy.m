%% load dataset
load 'data/gauss60_5+';
%% normalize dataset
mu = mean(dataset,2);
dataset = dataset - mu;
d_min = min(dataset,[],'all');
d_max = max(dataset,[],'all');
dataset = (dataset - d_min) ./ (d_max - d_min);
%% perform SVD
% [U,S,V] = svd(cov(dataset'));
[U,S,V] = svd(dataset);
diagS = diag(S);
Ntf = eigEnNum(diagS,0.9);figure;
semilogy(diagS(1:Ntf));

%% set gap parameters
% unmarred fraction
p = 0.5;
% unmarred number
n_unmarred = round(p * size(dataset,1));

%% set up masks
mask = makeMask(dataset,n_unmarred);
for j = 1:size(data_masked,2)
    M{j} = maskMatrix(mask(:,j));
end

%% marred average
data_masked = dataset .* mask;
mu_masked = sum(data_masked,2)./ (sum(mask,2) + eq(sum(mask,2),0));
data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* mu_masked;

max_iter = 10;
for num_iter = 1:max_iter
    [Um,Sm,~] = svd(data_masked);
    hold on;
    diagS = diag(Sm);
    Nt = eigEnNum(diagS,0.9);
%     Nt = 30;
    semilogy(diagS(1:Ntf));
    %% least squares problem
    A =  Um(:,1:Nt);
    a = zeros(Nt,size(data_masked,2));
    for j = 1:size(data_masked,2)
        B = data_masked(:,j);
        x = lscov(M{j}*A,M{j}*B);
        a(:,j) = x;
    end
    
    %%
    data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* (A*a);
end

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