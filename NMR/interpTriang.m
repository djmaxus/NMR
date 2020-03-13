function [Fq] = interpTriang(xq,Xv,F)
    %INTERPTRIANG Summary of this function goes here
    %   Detailed explanation goes here
    n = dsearchn(Xv,xq);
    gradF = calcGrad(Xv',F,'n',n);
    Fq = F(n,:) + (xq -Xv(n,:))*gradF;
end

% function [out] = calcGrad(X,F,varargin)
%     %% input processing
%     if length(varargin) == 2
%         if strcmp(varargin{1},'n')
%             n = varargin{2};
%         elseif strcmp(varargin{1},'x0')
%             x0 = varargin{2};
%             n = dsearchn(X',x0');
%         end
%     else
%         n = 1;
%     end
%     %% calculation
%     F = F - F(n,:);
%     X = X - X(:,n);
%     out = lscov(X',F);
% end



