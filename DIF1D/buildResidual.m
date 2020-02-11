function [Rt,Rx,Rxl,Rxr] = buildResidual(X,dt,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);
    %% time derivative
    Rt = (X(:,2)-X(:,1))./dt;
    
    %% right flux
    Rxr = zeros(N,1);
    if nmr
        mask = POD.Rxr.idx_mask;
    else
        mask = 1:N;
    end
    for i=1:length(mask)
        [flux]=calcFlux(X(:,2),D,bc,mask(i)+1);
        Rxr(mask(i)) = flux;
    end
    if nmr
        Rxr = Rxr - (POD.Rxr.mu + POD.Rxr.init);
        Rxr = POD.Rxr.O*Rxr;
        Rxr = forward(Rxr,POD.Rxr);
    end
    
    %% left flux
    Rxl = zeros(N,1);
    if nmr
        mask = POD.Rxl.idx_mask;
    else
        mask = 1:N;
    end
    for i=1:length(mask)
        [flux]=calcFlux(X(:,2),D,bc,mask(i));
        Rxl(mask(i)) = flux(1);
    end
    if nmr
        Rxl = Rxl - (POD.Rxl.mu + POD.Rxl.init);
        Rxl = POD.Rxl.O*Rxl;
        Rxl = forward(Rxl,POD.Rxl);
    end
    
    %% total flux
    Rx = Rxr - Rxl;
    
end

function [R] = calcFlux(X,D,bc,i)
    N = length(X);
    if bc.first
        if i==1
            delta = [X( i ),bc.val(1)];
        elseif i==N+1
            delta = [bc.val(2), X( i - 1 )];
        else
            delta = [X( i ),X( i - 1)];
        end
        x = max(delta(2),delta(1));
        R = -D( i ) * (x^2+1) * ( delta(1) - delta(2) );
    else
        if i==1
            R = bc.val(1);
        elseif i==N+1
            R = bc.val(2);
        else
            x = max(X(i),X(i-1));
            R = -D( i ) * (x^2+1) * ( X( i ) - X( i - 1 ) );
        end
    end
end
