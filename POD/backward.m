function [Y] = backward(X,POD)
    %FORWARD Summary of this function goes here
    %   Detailed explanation goes here
    Y = POD.Ut'*(X - POD.init - POD.mu);
end

