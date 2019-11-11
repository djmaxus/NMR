function [POD] = makePOD(data, frac_energy)
    %MAKEPOD Summary of this function goes here
    %   Detailed explanation goes here
    [POD.U,POD.S,POD.V] = svd(data);
    POD.N = eigEnNum(diag(POD.S),frac_energy);
    POD.Ut = POD.U(:,1:POD.N);
    
    f = figure;
    f.WindowState = 'maximized';
    en_frac = diag(POD.S)./sum(diag(POD.S));
    
    semilogy(en_frac);
    
    hold on;
    xline(POD.N);
    
    title('POD');
    legend('eigen value',...
           ['reduction, frac = ',num2str(frac_energy)]);
    xlabel('Number of POD components');
    waitforbuttonpress;
end

