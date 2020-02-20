function [POD] = doNMR(FOM,fX,fR,fS)
    %DONMR Summary of this function goes here
    %   Detailed explanation goes here
    
    %% X POD
    POD.X = makePOD(FOM.X,true,true,fX);
    %% Residual POD
    POD.Rx = makePOD(FOM.Rx,true,true,fR);
    POD.Rx = makeOblique(fS,POD.Rx);
    %% Jacobian
    POD.Jnz = makePOD(FOM.Jnz,false,true,fR);
    POD.Jnz = makeOblique(fS,POD.Jnz);
    
    f = figure;
    f.WindowState = 'maximized';
    Nd = 0;
    names = fieldnames(POD);
    for basis = {POD.X,POD.Rx,POD.Jnz}
        F = basis{1};
        en_frac = diag(F.S)./sum(diag(F.S));
        
        semilogy(en_frac);
        Nd = max(Nd,min(size(F.S)));
        xlim([1 Nd]);
        
        hold on;
        plot(F.N,en_frac(F.N),'o');
          
        hold on;
    end
    xlabel('Number of POD components');
    ylabel('Component energy');
    legend('U',['N_{',names{1},'}'],'Residual',['N_{',names{2},'}',],'Jacobian',['N_{',names{3},'}']);
    title('POD');
    hold off;
end
