function [out] = genRom1(param)
    %TRAINROM1 Summary of this function goes here
    %   Detailed explanation goes here
    %%
    Ns = size(param.vals,1);
    sol = cell(Ns,1);
    parfor i = 1:Ns
        param_i = param;
        param_i.bc.val = param.vals(i,:);
        sol{i} = solveDif(param_i.Nx,param_i.Nt,param_i,'train_mode',true);
        disp(i);
    end
    out = sol;
    %%
    %     if length(varargin) >=1
    %         out = varargin{1};
    %     else
    %         out = cell(size(param.grids));
    %         for iNx = 1:length(param.grids)
    %             out{iNx}.Nx = 0;
    %         end
    %     end
    %
    
    %
    %     for iNx = 1:length(param.grids)
    %         found_grid = false;
    %         for i = 1:length(out)
    %             if out{iNx}.Nx == param.grids(iNx)
    %                 found_grid = true;
    %                 break;
    %
    %             end
    %         end
    %         if found_grid
    %             continue;
    %         else
    %             data.Nx = param.grids(iNx);
    %             disp(['Nx: ',num2str(data.Nx)]);
    %
    %             out{iNx} = data;
    %         end
    %     end
end

