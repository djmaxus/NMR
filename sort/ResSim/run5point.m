function [wellState,resState] = run5point(perm)
%RUN5POINT Summary of this function goes here
%   Detailed explanation goes here
problem = fivePoint(perm);
[wellState, resState] = ResSimSPCA(problem);
end

