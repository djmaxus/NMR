function [] = plotMSE(data,POD,mask,grd)
    %PLOTMSE Summary of this function goes here
    %   Detailed explanation goes here
    out.train = calcMSE(POD.U,POD.S,data.train);
    out.test = calcMSE(POD.U,POD.S,data.test);
    
    delta = data.test - data.restored.test;
    out.restored = sum(delta .^ 2,'all')/numel(data.test);
    
    delta = data.test - data.res_grd.test;
    out.grd = sum(delta .^ 2,'all')/numel(data.test);
    
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
    
    title('Mean Square Error');
    legend('train','test',['restored, p = ',num2str(mask.frac)],...
                          ['greedy p = ',num2str(grd.frac)]);
    xlabel('Number of POD components');
end

