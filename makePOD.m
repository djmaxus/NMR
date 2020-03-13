function [POD,fig] = makePOD(data, init, mu, param,varargin)
    %MAKEPOD Summary of this function goes here
    %   Detailed explanation goes here
    if length(varargin) >= 1
        viz = true;
        fig = varargin{1};
        fig.Name = 'PCA';
    else
        viz = false;
    end
    %% Substract initial
    if init
        POD.init = data(:,1);
        data = data - POD.init;
    else
        POD.init = zeros(size(data(:,1)));
    end
    
    
    %% Substract mean
    if mu
        POD.mu = mean(data,2);
        data = data - POD.mu;
    else
        POD.mu = zeros(size(data(:,1)));
    end
    
    %% Crop Data
    %     cropData;
    %% SVD
    [POD.U,S,~] = svd(data,'econ');
    Nd = size(S,1);
    dS = diag(S)/sum(diag(S));
    if param > 1
        POD.N = min(round(param),Nd);
        
        POD.en = sum(dS(1:POD.N));
    elseif param > 0
        POD.en = param;
        POD.N = eigEnNum(dS,param);
    else
        error('incorrect parameter');
    end
    %% visualize
    if viz
        % principal components
        subplot(1,2,1);
        for i = 1:size(data,2)
            plot(data(:,i),'k-');
            hold on;
        end
        % principal components
        xmax = max(abs(data),[],'all');
        idx = [1,2,Nd];
        stl = {'r-','b-','g-'};
        for i = 1:length(idx)
            waitforbuttonpress;
            u = POD.U(:,idx(i));
            umax = max(abs(u),[],'all');
            u = u .* xmax ./ umax ./ i;
            plot(u,stl{i},'LineWidth',1);
            hold on;
        end
        hold off;
        % reconstruction error
        X = POD.U * POD.U' * data;
        err = norm(data - X);
        title(['Reconstruction error: ',num2str(err)]);
        drawnow;
        % energy
        waitforbuttonpress;
        subplot(1,2,2);
        enS = dS;
        for i = 1:Nd
            enS(i) = sum(dS(1:i));
        end
        semilogx(1:Nd,enS);
        xline(POD.N);
        yline(POD.en);
        ylim([0,1]);
        title([num2str(100 * POD.en),' % of basis energy',' in ',num2str(POD.N),' components']);
        drawnow;
    end
    
end
