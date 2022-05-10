function [F] = genField(N,Fmin,Fmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
F = rsgengXY(1,N,1,1,0.2);
F = Fmin + (Fmax - Fmin) * (F - min(F)) / (max(F) - min(F));
end

