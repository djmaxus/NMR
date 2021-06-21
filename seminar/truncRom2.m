function [POD,R,J] = truncRom2(POD,train)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train)
        for j = 1:length(train{i}.nlin)
            for m = 1:(train{i}.nlin{j}.iters + 1)
                R = [R,train{i}.nlin{j}.residual{m}.W];
                J = [J,train{i}.nlin{j}.jacobian{m}.Jf.*[1,-1]];
            end
        end
    end
    
    r0 = zeros(size(R,1),1);
    
    POD.Rx = PCA(R,0.99,'init',r0,'sval_deg',1,'robust',false);
    POD.Jnz = PCA(J,0.99,'sval_deg',1);
    
    Nr = max(POD.Rx.npc(),POD.X.npc());
    Nj = max(POD.Jnz.npc(),POD.X.npc());
    
    POD.Rx = POD.Rx.setNPC(Nr);
    POD.Jnz = POD.Jnz.setNPC(Nj);
end

