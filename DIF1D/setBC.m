function [bc] = setBC(kind,pMax)
    %SETBC Summary of this function goes here
    %   Detailed explanation goes here
    bc.first = -kind + 2;
    bc.val = rand(1,2) * 2 * pMax - pMax;
end

