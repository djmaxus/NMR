function [A] = applyOblique(B,POD)
    %RECONSTRUCT Summary of this function goes here
    %   Detailed explanation goes here
    %   Applies oblique projection assuming there are no initial nor mean value calculated for POD
    A = POD.O * B;
%     disp(sum(POD.mask'*abs(A-B),'all'));
%     if size(B,1) == size(B,2)
%         Dc = diag(A-B);
%         disp(sum(abs(Dc),'all'));
%         Dd = diag(A-B,-1);
%         disp(sum(abs(Dd),'all'));
%         Du = diag(A-B,1);
%         disp(sum(abs(Du),'all'));
%     end
%     disp(sum(abs(A-B),'all'));
end
