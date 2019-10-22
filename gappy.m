%% load dataset
load 'data/gauss60_5+';
%% normalize dataset
mu = mean(dataset,2);
dataset = dataset - mu;
%% perform SVD
[U,S,V] = svd(dataset);
%% set gap parameters
% unmarred fraction
p = 0.3;
% unmarred number
n_unmarred = round(p * size(dataset,1));
%% set up masks
i_unmarred = randi([1,size(dataset,1)],n_unmarred,size(dataset,2));
mask = full(makeMask(dataset,n_unmarred));
data_masked = dataset .* mask;
%% marred average
diagS = diag(S);
Nt = eigEnNum(diagS,0.9);
semilogy(diagS(1:Nt));
mu_masked = sum(data_masked,2)./ (sum(mask,2) + eq(sum(mask,2),0));
data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* mu_masked;
max_iter = 20;
for num_iter = 1:max_iter
    [Um,Sm,~] = svd(data_masked);
    hold on;
    diagS = diag(Sm);
    %     Nt = eigEnNum(diagS,0.8);
    semilogy(diagS);
    %% least squares problem
    A =  Um(:,1:Nt);
    a = A'*data_masked;
    for j = 1:size(data_masked,2)
        M = maskMatrix(mask(:,j));
        B = data_masked(:,j);
        x = lscov(M*A,M*B);
        delta = norm(x - a(:,j));
        a(:,j) = x;
    end
    
    %%
    data_masked = data_masked .* eq(mask,1) + eq(mask,0) .* (A*a);
end
figure; 
histogram(dataset(:,1),50);
hold on;
histogram(data_masked(:,1),50);
%% check if entry is non-unique in array
function [flag] = is_repeat(entry_index,array)
sum_left = sum(eq(array(entry_index),array(1:entry_index-1)));
sum_right = sum(eq(array(entry_index),array(entry_index+1:end)));
flag = (sum_left + sum_right) > 0;
end

%% calculate mask
function mask = makeMask(dataset,n_unmarred)
mask = spalloc(size(dataset,1),size(dataset,2),size(dataset,2)*n_unmarred);
i_unmarred = randi([1,size(dataset,1)],n_unmarred,size(dataset,2));
for i = 1:size(dataset,2)
    for j = 2:n_unmarred
        while is_repeat(j,i_unmarred(:,i))
            val = randi([1,size(dataset,1)],1,1);
            i_unmarred(j,i) = val;
        end
        mask(i_unmarred(j,i),i) = 1;
    end
end
end
function out = maskMatrix(mask)
out = spalloc(sum(mask,1),size(mask,1),sum(mask,1));
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