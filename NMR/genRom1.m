function [out] = genRom1(param,varargin)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >=1
        out = varargin{1};
    else
        out = cell(size(param.grids));
        for iNx = 1:length(param.grids)
            out{iNx}.Nx = 0;
        end
    end
    
    Ns = size(param.vals,1);
    sol = cell(Ns,1);
    
    for iNx = 1:length(param.grids)
        found_grid = false;
        for i = 1:length(out)
            if out{iNx}.Nx == param.grids(iNx)
                found_grid = true;
                break;
                
            end
        end
        if found_grid
            continue;
        else
            data.Nx = param.grids(iNx);
            D = setDiff(data.Nx,param.D);
            disp(['Nx: ',num2str(data.Nx)]);
            parfor i = 1:Ns
                bc_i = param.bc;
                bc_i.val = param.vals(i,:);
                sol{i} = solveDif(data.Nx,param.Nt,bc_i,D);
                disp(i);
            end
            data.sol = sol;
            out{iNx} = data;
        end
    end
end

