function [POD] = makeOblique(frac,POD)
    %MAKEOBLIQUE Summary of this function goes here
    %   Detailed explanation goes here
    %% greedy mask
    Nx = min(size(POD.U));
    Nm = round(frac*Nx);
    Nm = max(Nm,POD.N);
    mask = greedy(POD.U(:,1:Nm));
    %% oblique projector
    idx_mask = zeros(Nm,1);
    for i = 1:Nm
        for j=1:Nx
            if mask(j,i) ==1
                idx_mask(i) = j;
                break;
            end
        end
    end
    %% add boundary cells
%     find_left = min(idx_mask) <= 1;
%     find_right = max(idx_mask) >= Nx;
%     if ~find_left
%         idx_mask = [idx_mask;1];
%         dmask = zeros(Nx,1);
%         dmask(1) = 1;
%         mask = [mask,dmask];
%     end
%     if ~find_right
%         idx_mask = [idx_mask;Nx];
%         dmask = zeros(Nx,1);
%         dmask(Nx) = 1;
%         mask = [mask,dmask];
%     end
    %% set data
    POD.mask = mask;
    POD.idx_mask = idx_mask;
    %% calculate oblique projector
    POD.O = oblique(mask,POD.U(:,1:POD.N));
end
