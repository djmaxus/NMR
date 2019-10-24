function out = eigEnNum(S,fraction)
total = sum(S);
for i = 1:length(S)
    if sum(S(1:i)) >= fraction * total
        out = i;
        break;
    end
end
end