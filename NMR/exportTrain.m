function [] = exportTrain(data,path)
    %EXPORTTRAIN Summary of this function goes here
    %   Detailed explanation goes here
    formatSpec = '%e;';
    Ng = length(data);
    for i = 1:Ng
        Nx = data{i}.Nx;
        fname = [path,'dataset_',num2str(Nx),'.csv'];
        fid = fopen(fname,'w');
        if fid == -1
            error(['Cannot open file ',fname]);
        else
            Ns = length(data{i}.sol);
            for j = 1:Ns
                Nt = size(data{i}.sol{j}.X,2);
                for k = 1:Nt
                    for  m = 1:Nx
                        fprintf(fid,formatSpec,data{i}.sol{j}.X(m,k));
                    end
                    fprintf(fid, '\n');
                end
            end
        end
        fclose(fid);
    end
end

