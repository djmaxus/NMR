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

function [flag] = is_repeat(entry_index,array)
    sum_left = sum(eq(array(entry_index),array(1:entry_index-1)));
    sum_right = sum(eq(array(entry_index),array(entry_index+1:end)));
    flag = (sum_left + sum_right) > 0;
end