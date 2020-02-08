function [Rt,Rx,Jt,Jx] = buildLinsystem(X,dt,D,bc,mask)
    N = size(X,1);
    Jt = zeros(N,1);
    % time derivative
    Rt = (X(:,2) - X(:,1))./dt;
    for i=1:N
        Jt(i,i) = 1/dt;
    end

    % flux
    Rx = zeros(N,1);
    %TODO: make 3 diags
    Jx = spalloc(N,N,3*N-2);
    for i=mask
        flux = [0;0];
        Jf = [0,0;0,0];
        [flux(1),Jf(1,:)]=calcFlux(X,D,bc,i);
        [flux(2),Jf(2,:)]=calcFlux(X,D,bc,i+1);

        Rx(i) = flux(2) - flux(1);
        Jx(i,i-1) = - Jf(1,1);
        Jx(i,i) = Jf(2,1)-Jf(1,2);
        Jx(i,i+1) = Jf(2,2);
    end
end

function [R,J] = calcFlux(X,D,bc,i)
    if i==1
        R = bc(1);
    elseif i==N+1
        R = bc(2);
    else
        x = max(X(i),X(i-1));
        imax = eq(x,X(i-1)) + eq(x,X(i))*2;
        R = D( i ) * x^2 * ( X( i ) - X( i - 1 ) );
    end

    J =[0,0];
    for j=1:2
        if and(i==1,j==1)
            J( j ) = 0;
        elseif and(i==N+1,j==2)
            J( j ) = 0;
        else
            J( j ) = D( j ) * (2 * x * eq( imax, j ) * (X( j ) - X( j - 1 )) + x ^ 2 * (eq( j, 2 ) - eq( j, 1 )));
        end
    end
end
