function [POD] = truncRom2(POD,train)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train)
        R = [R,train{i}.Rx];
        J = [J,train{i}.Jnz];
    end
%     R = train.Rx;
%     J = train.Jnz;
    r0 = zeros(size(R,1),1);
    j0 = zeros(size(J,1),1);
    POD.Rx = PCA(R,0.999,'init',r0,'sval_deg',1,'robust',true);
    POD.Jnz = PCA(J,0.999,'sval_deg',1);
    
    Nr = max(POD.Rx.npc(),POD.X.npc());
    Nj = max(POD.Jnz.npc(),POD.X.npc());
    
    POD.Rx = POD.Rx.setNPC(Nr);
    POD.Jnz = POD.Jnz.setNPC(Nj);
end

