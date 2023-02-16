classdef BuilderPOD < handle
    
    properties
        data;
        shift;
        
        u;
        s;
        v;
        
        w;
        
        rk (1,1) = 0;
        
        gappy_mask (1,:) {mustBeInteger, mustBePositive};
        gappy_projection (:,:) {sparse};
        
        mask_extension(1,:) {mustBeInteger, mustBePositive} = [];
    end
    
    methods
        function obj = BuilderPOD(data,shift)
            obj.data = data;
            obj.shift = shift;
            
            [obj.u, s, obj.v] = svd(data - shift,'econ');
            
            variance = sqrt(var(data,1,'all'));
            obj.w = norm(shift) + abs(variance);
            
            s = diag(s);
            s = s./(size(data,2)-1);
            
            obj.s = s;
        end
        
        function [obj] = set_rank(obj,rk)
            rk = ceil(rk);
            max_rank = size(obj.u, 2);
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
        
        function [obj] = set_lower_bound(obj,lower_bound)
            obj.set_rank(max([obj.rk,lower_bound]));
        end
        
        function [encoder, obj] = build_encoder(obj)
            encoder = EncoderUnitary(obj.u(:,1:obj.rk),obj.shift);
        end
        
        function [encoder, obj] = build_encoder_scaled(obj)
            u_scaled = obj.u(:,1:obj.rk).*obj.w;
            encoder = EncoderLinear(u_scaled,obj.shift);
        end
        
        function [obj] = set_gappy_mask(obj,increase_factor)
            arguments
                obj (1,1) BuilderPOD
                increase_factor (1,1) {mustBeReal,mustBePositive} = 1;
            end
            increase_factor = max(increase_factor,1);
            mask_size = ceil(increase_factor * obj.rk);
            [obj.gappy_mask] = calc_gappy_mask(obj.u,mask_size);
            obj.mask_extension = [];
        end
        
        function [obj] = set_mask_extension(obj, indices)
            arguments
                obj (1,1) BuilderPOD
                indices (1,:) {mustBeInteger, mustBePositive}
            end
            
            max_index = size(obj.data,1);
            indices = min(indices,max_index);
            indices = max(indices,1);
            obj.mask_extension = unique([obj.mask_extension, indices]);
        end
        
        function [obj] = set_mask_extension_normalized(obj, range)
            arguments
                obj (1,1) BuilderPOD
                range (1,:) double
            end
            
            a = (size(obj.data,1) - 1);
            range = a .* range + 0.5; % make it from [0.5,N-0.5]
            range = ceil(range); % shift it from [1,N]
            obj.mask_extension = [obj.mask_extension, range];
        end
        
        function [encoder_gappy, obj] = build_gappy(obj)
            encoder_pod = obj.build_encoder();
            gappy_mask_compound = unique([obj.gappy_mask,obj.mask_extension]);
            [~,gappy_decode_matrix] = oblique_projector(gappy_mask_compound,encoder_pod.decoding_transform);
            encoder_gappy = EncoderGappy(gappy_decode_matrix,encoder_pod.shift,gappy_mask_compound);
        end
    end
end
