function [A] = applyOblique(B,POD)
    %RECONSTRUCT Summary of this function goes here
    %   Detailed explanation goes here
    %   Applies oblique projection assuming there are no initial nor mean value calculated for POD
    A = POD.O * (B-POD.init - POD.mu)+POD.init + POD.mu;
    disp(sum(POD.mask'*abs(A-B),'all'));
%     disp(sum(abs(A-B),'all'));
end
