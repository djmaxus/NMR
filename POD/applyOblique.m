function [A] = applyOblique(B,POD)
    %RECONSTRUCT Summary of this function goes here
    %   Detailed explanation goes here
    
    A = POD.O * B;
end

