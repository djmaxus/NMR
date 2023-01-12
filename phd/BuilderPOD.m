classdef BuilderPOD
    
    properties
        data;
        shift;
        
        u;
        s;
        v;
        
        rk;
    end
    
    methods
        function obj = BuilderPOD(data,shift)
            obj.data = data;
            obj.shift = shift;
            
            [obj.u, s, obj.v] = svd(data - shift,'econ');
            
            obj.s = diag(s);   
        end
        
        function [obj] = set_rank(obj,rk)
            rk = ceil(rk);
            max_rank = size(obj.data, 2);
            obj.rk = max(min(rk,max_rank),1);
        end
        
        function [obj] = set_rank_fraction(obj,fraction)
            fraction = min(max(fraction,0),1);
            [obj] = obj.set_rank(size(obj.data,2)*fraction);
        end
        
        function [obj] = set_energy_residual(obj,energy_residual)
            partial_residuals = cumsum(obj.s(end:-1:1));
            total_energy = partial_residuals(end);
            
            scaled_residual = energy_residual * total_energy;
            energy_rank = length(partial_residuals(partial_residuals > scaled_residual))+1;
            
            obj = obj.set_rank(energy_rank);
        end
        
        function obj = set_energy(obj,energy)
            obj = obj.set_energy_residual(1 - energy);
        end
        
        function [obj,encoder] = build_encoder(obj)
            encoder = EncoderUnitary(obj.u(:,1:obj.rk),obj.shift);
        end
    end
end

