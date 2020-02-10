function [POD,mask,idx_mask] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    Nx = size(FOM.X,1);
    %% X POD
    POD.X = makePOD(FOM.X,fX);
    %% Residual POD
    POD.Rx = makePOD(FOM.Rx,fR);
    %% greedy mask
    frac = fS;
    Nm = round(frac*Nx);
    Nm = max(Nm,POD.Rx.N);
    mask = greedy(POD.Rx.U(:,1:Nm));
    % mask = mask(:,1:Nm);
    %% oblique projector
    POD.Rx.O = oblique(mask,POD.Rx.Ut);
    idx_mask = zeros(Nm,1);
    for i = 1:Nm
        for j=1:Nx
            if mask(j,i) ==1
                idx_mask(i) = j;
                break;
            end
        end
    end
    find_left = false;
    find_right = false;
    for i = idx_mask'
        if i == 1
            find_left = true;
        elseif i == Nx
            find_right = true;
        end
    end
    if ~find_left
        idx_mask = [idx_mask;1];
        dmask = zeros(Nx,1);
        dmask(1) = 1;
        mask = [mask,dmask];
    end
    if ~find_right
        idx_mask = [idx_mask;Nx];
        dmask = zeros(Nx,1);
        dmask(Nx) = 1;
        mask = [mask,dmask];
    end
end

