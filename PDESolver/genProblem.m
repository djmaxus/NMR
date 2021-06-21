function [problem] = genProblem(name,Nx,param,x0,varargin)
    %GENPROBLEM Summary of this function goes here
    %   Detailed explanation goes here
    
    % encoder
    encode = @(x) x;
    decode = @(x) x;
    if ~isempty(varargin)
        encode = varargin{1}.encode;
        decode = varargin{1}.decode;
        varargin(1) = [];
    end
    
    % equations
    res = @(x,dt) buildResidual(decode(x),dt,param,varargin{:});
    if ~isempty(varargin)
        varargin(1) = [];
    end
    
    jac = @(x,dt) buildJacobian(decode(x),dt,param,varargin{:});
    
    % y0
    y0 = encode(x0);
    
    % data
    problem = struct(...
        'name'    , name,...
        'Nx'      , Nx,...
        'residual', res,...
        'jacobian', jac,...
        'x0'      , y0,...
        'param'   , param...
        );
end

