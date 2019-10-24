function [] = plotMSE(data,POD)
    %PLOTMSE Summary of this function goes here
    %   Detailed explanation goes here
    out.train = calcMSE(POD.U,POD.S,data.train);
    out.test = calcMSE(POD.U,POD.S,data.test);
    f = figure;
    f.WindowState = 'maximized';
    semilogy(out.train);
    hold on;
    semilogy(out.test);
    hold on;
    xline(POD.N);
    title('Mean Square Error');
    legend('train','test');
    xlabel('Number of POD components');
end

