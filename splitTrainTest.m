function [data] = splitTrainTest(source,frac_test)
    %SPLITTRAINTEST Summary of this function goes here
    %   Detailed explanation goes here
    N = round(size(source,2)*frac_test);
    i_test = randi([1,size(source,2)],N,1);
    i_test = sort(unique(i_test));
    test = zeros(size(source,1),length(i_test));
    train = source;
    for i = 1:length(i_test)
        test(:,i) = source(:,i_test(i));
        train(:,i_test(i)-i) = [];
    end
    data.train = train;
    data.test = test;
end

