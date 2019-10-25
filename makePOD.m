function [POD] = makePOD(data, frac_energy)
    %MAKEPOD Summary of this function goes here
    %   Detailed explanation goes here
    [POD.U,POD.S,POD.V] = svd(data);
    POD.N = eigEnNum(diag(POD.S),frac_energy);
    POD.Ut = POD.U(:,1:POD.N);
end

