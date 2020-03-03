function [Y] = backward(X,POD)
    %FORWARD Summary of this function goes here
    %   Detailed explanation goes here
%     Y = (X - POD.mu - POD.init)./POD.sigma;
    Y = POD.U(:,1:POD.N)'*(X - POD.init - POD.mu);
end
