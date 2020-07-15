function [POD] = truncRom2(POD,train)
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
    POD.Jnz = PCA(J,0.99);
    
    Nr = max(POD.Rx.npc(),POD.X.npc());
    Nj = max(POD.Jnz.npc(),POD.X.npc());
    
    POD.Rx = POD.Rx.setNPC(Nr);
    POD.Jnz = POD.Jnz.setNPC(Nj);
end

