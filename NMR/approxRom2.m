function [POD] = approxRom2(POD,train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train.sol)
        R = [R,train.sol{i}.Rx];
        J = [J,train.sol{i}.Jnz];
    end
    
    add_frac = 0.5;
    %% Res POD + greedy
    nR = round(POD.Rx.npc() * (1 + add_frac));
    [~,mask_R] = greedy(POD.Rx.u(:,1:POD.Rx.npc()),nR);
    %% Jac POD + greedy
    nJ = round(POD.Jnz.npc() * (1 + add_frac));
    [~,mask_J] = greedy(POD.Jnz.u(:,1:POD.Jnz.npc()),nJ);
    %% Common mask
    [idx,mask] = mergeMasks(mask_R,mask_J);
    POD.Rx.idx_mask = idx;
    POD.Rx.mask = mask;
    POD.Jnz.idx_mask = idx;
    POD.Jnz.mask = mask;
    %%
%     POD.Rx = POD.Rx.setNPC(length(idx));
%     POD.Jnz = POD.Jnz.setNPC(length(idx));
    %% Res oblique
    POD.Rx.O = oblique(mask,POD.Rx.u(:,1:POD.Rx.npc()));
    [~,dR] = applyOblique(R,POD.Rx);
    disp(dR);
    %% Jac oblique
    POD.Jnz.O = oblique(mask,POD.Jnz.u(:,1:POD.Jnz.npc()));
    [~,dJ] = applyOblique(J,POD.Jnz);
    disp(dJ);
end

function [idx_mask,P] = mergeMasks(mask_1,mask_2)
    mask = sum(mask_1,2) + sum(mask_2,2);
    M = sum(ne(mask,0),'all');
    idx_mask = zeros(1,M);
    P = spalloc(size(mask,1),M,M);
    idx = 1;
    for i = 1:length(mask)
        if mask(i) ~= 0
            idx_mask(idx) = i;
            P(i,idx) = 1;
            idx = idx+1;
        end
    end
end
