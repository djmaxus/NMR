function [] = plotMSE(data,POD,mask,grd)
    %PLOTMSE Summary of this function goes here
    %   Detailed explanation goes here
    out.train = calcMSE(POD.U,POD.S,data.train);
    out.test = calcMSE(POD.U,POD.S,data.test);
    
    delta = data.test - data.restored.test;
    delta = diag(delta'*delta);
    delta = sum(delta)/length(delta);
    out.restored = delta;
    
    delta = data.test - data.res_grd.test;
    delta = diag(delta'*delta);
    delta = sum(delta)/length(delta);
    out.grd = delta;
    
    f = figure;
    f.WindowState = 'maximized';
    semilogy(out.train);
    
    hold on;
    semilogy(out.test);
    
    hold on;
    plot(POD.N,out.restored,'o');
    
    hold on;
    plot(POD.N,out.grd,'o');
    
    hold on;
    xline(POD.N);
    
    title('Mean Squared Error');
    legend('train','test',['restored, p = ',num2str(mask.frac)],...
                          ['greedy p = ',num2str(grd.frac)]);
    xlabel('Number of POD components');
end

