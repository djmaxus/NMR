function [Rt,W] = buildResidual(X,dt,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);
    %% time derivative
    Rt = (X(:,2)-X(:,1))./dt;

    %% fluxes
    
    mask = 1:N+1;
    if nmr
        mask = POD.idx_mask;
    end
    W = zeros(N+1,1);
    for i = mask
        W(i) = calcFlux(X(:,2),D,bc,i);
    end
    if nmr
        W = applyOblique(W,POD);
    end 
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
        mult = D.model(x);%(x^2+1);
        R = -D.val( i ) * mult(1) * ( delta(1) - delta(2) );
    else
        if i==1
            R = bc.val(1);
        elseif i==N+1
            R = bc.val(2);
        else
            x = max(X(i),X(i-1));
            R = -D.val( i ) * (x^2+1) * ( X( i ) - X( i - 1 ) );
        end
    end
end
