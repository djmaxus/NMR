function [POD] = approxRom2(POD,train,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    R = [];
    J = [];
    for i = 1:length(train)
        R = [R,train{i}.Rx];
        J = [J,train{i}.Jnz];
    end
%     add_frac = 1.0;
    %% Res POD + greedy
    nR = round(size(R,1)*0.01);
    nR = max(nR,POD.Rx.npc());
%     nR = round(POD.Rx.npc()*(1+add_frac));
    [~,mask_R] = greedy(POD.Rx.u(:,1:POD.Rx.npc()),nR);
    
    %% Jac POD + greedy
    nJ = round(size(J,1)*0.01);
    nJ = max(nJ,POD.Jnz.npc());
%     nJ = round(POD.Jnz.npc()*(1+add_frac));
    [~,mask_J] = greedy(POD.Jnz.u(:,1:POD.Jnz.npc()),nJ);
    
    %% Common mask
    [idx,mask] = mergeMasks(mask_R,mask_J);
    
    POD.Rx.mask = mask;
    POD.Jnz.mask = mask;
    
    POD.Rx.idx_mask = idx;
    POD.Jnz.idx_mask = idx;
    
    %% Res oblique
    uR = POD.Rx.u(:,1:POD.Rx.npc());
    uJ = POD.Jnz.u(:,1:POD.Jnz.npc());
    POD.Rx.O = oblique(mask,uR);
    POD.Jnz.O = oblique(mask,uJ);
    
    %% display Gappy scores
    uR = mask' * uR;
    uJ = mask' * uJ;
    MR = uR' * uR;
    MJ = uJ' * uJ;
    POD.score.R = rcond(MR);
    POD.score.J = rcond(MJ);
    
    Rg = applyOblique(R,POD.Rx);
    Jg = applyOblique(J,POD.Jnz);
    dR = sqrt(meansqr(Rg-R));
    dJ = sqrt(meansqr(Jg-J));
    disp(...
        [...
        '___Gappy scores___',newline(),...
        '  Res: ', num2str(dR,'%6.2e'),newline(),...
        '  Jac: ', num2str(dJ,'%6.2e'),newline(),...
        '------------------'...
        ]...
        );
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
