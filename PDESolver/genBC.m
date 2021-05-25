function [bc] = genBC(kind,lim)
    %BC Summary of this function goes here
    %   Detailed explanation goes here
    bc.kind = kind;
    bc.val = rand(1,2) * (lim(2)-lim(1)) + lim(1);
end

