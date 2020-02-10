function [Rt,Rx] = buildResidual(X,dt,D,bc,mask,POD)
    N = size(X,1);
    nmr = length(mask) < N;
    % time derivative
    Rt = (X(:,2)-X(:,1))./dt;

    % flux
    Rx = zeros(N,1);
    for i=1:length(mask)
        flux = [0;0];
        [flux(1)]=calcFlux(X(:,2),D,bc,mask(i));
        [flux(2)]=calcFlux(X(:,2),D,bc,mask(i)+1);

        Rx(mask(i)) = flux(2) - flux(1);
    end
    if nmr
        Rx = Rx - (POD.Rx.mu + POD.Rx.init);
        Rx = POD.Rx.O*Rx;
        Rx = forward(Rx,POD.Rx);
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
