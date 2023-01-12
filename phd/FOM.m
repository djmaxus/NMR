classdef FOM
    %FOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dx;
        ul;
        ur;
        u0;
        kappa;
        q;
        calc_g = @(u) [1 + 0.5 * u.^2, u];
    end
    
    properties(Constant)
        
    end
    
    methods
        function obj = FOM(Nu, ul, ur, gen_initial, gen_kappa, source)
            Nu = max([Nu,3]);
            Nf = Nu - 1;
            obj.dx = Nf^(-1);
            
            obj.ul = ul;
            obj.ur = ur;
            
            obj.u0 = gen_initial(Nu);
            obj.u0(1) = ul;
            obj.u0(end) = ur;
            
            obj.kappa = gen_kappa(Nf);
            
            pos = max(min(round( (Nu-1)*source.pos + 1 ), Nu-1), 2);
            obj.q = sparse(pos,1,source.value / obj.dx, Nu, 1);
        end
        
        function [residual, jacobian] = calc_linearization(obj,u0,dt,u)
            [residual, jacobian] = pde_linearization(u0, dt, u, obj.dx, obj.ul, obj.ur, obj.kappa, obj.q, obj.calc_g);
        end
        
        function plot_problem(obj,ax)
            todo(); % plotProblem(axes(figure),FOM_test);
        end
        
        function plot_solution(obj,where)
            todo(); % plotSol(figure(),FOM_test);
        end
    end
end

