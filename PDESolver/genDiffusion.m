function [D] = genDiffusion(Nx,amp,per,scale,varargin)
    %SETDIFF Summary of this function goes here
    %   Detailed explanation goes here
    
    D.model = @(u)([1 + u.^2, 2*u]);
    
    if length(varargin) == 1
        D.model = varargin{1};
    end
    
    X = (0:Nx)'./ (Nx+1);
    
    circ = 2 * pi;
    
    D.val = scale * (amp * sin(per * circ * X) + 1);
end

