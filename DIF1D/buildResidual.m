function [Rt,Rx] = buildResidual(X,dt,D,bc,mask)
    N = size(X,1);
    % time derivative
    Rt = (X(:,2) - X(:,1))./dt;

    % flux
    Rx = zeros(N,1);
    for i=1:length(mask)
        flux = [0;0];
        [flux(1)]=calcFlux(X(:,2),D,bc,mask(i));
        [flux(2)]=calcFlux(X(:,2),D,bc,mask(i)+1);

        Rx(i) = flux(2) - flux(1);
    end
end

function [R] = calcFlux(X,D,bc,i)
    N = length(X);
    if i==1
        R = bc(1);
    elseif i==N+1
        R = bc(2);
    else
        x = max(X(i),X(i-1));
        R = D( i ) * x^2 * ( X( i ) - X( i - 1 ) );
    end
end
