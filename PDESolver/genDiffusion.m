function [D] = genDiffusion(Nx,amp,per,scale,varargin)
    %SETDIFF Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) == 1
        D.model = model;
    else
        D.model = @(u)([max(u.^2,0)+1,max(2*u,0)]);
    end
    D.val = scale * (amp * sin(per*2*pi*(1:Nx+1)'./(Nx+1)) + 1);
end

