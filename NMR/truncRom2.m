function [POD] = truncRom2(POD,train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train.sol)
        R = [R,train.sol{i}.Rx];
        J = [J,train.sol{i}.Jnz];
    end
    r0 = zeros(size(R,1),1);
    j0 = zeros(size(J,1),1);
    POD.Rx = PCA(R,0.99,'init',r0);
    POD.Jnz = PCA(J,0.99,'init',j0);
    
    Nr = max(POD.Rx.npc(),POD.X.npc());
    Nj = max(POD.Jnz.npc(),POD.X.npc());
    
    POD.Rx.crit_num = true;
    POD.Rx.param = Nr;
    POD.Jnz.crit_num = true;
    POD.Jnz.param = Nj;
end

