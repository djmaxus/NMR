function [X] = forward(Y,POD)
    %FORWARD Summary of this function goes here
    %   Detailed explanation goes here
    X =  POD.init + POD.mu + POD.U(:,1:POD.N)*Y;
end
