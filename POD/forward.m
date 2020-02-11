function [X] = forward(Y,POD)
    %FORWARD Summary of this function goes here
    %   Detailed explanation goes here
    X =  POD.Ut*Y;
end

