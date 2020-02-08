function [Jt,Jx] = buildJacobian(X,dt,D,mask)
    N = size(X,1);
    % time derivative
    Jt = ones(N,1)./dt;

    % flux
    Jx = spalloc(N,N,3*N-2);
    %TODO: make 3 diags
    Ju = zeros(N-1);
    Jc = zeros(N);
    Jd = zeros(N-1);
    for i=mask
        Jf = [0,0;0,0];
        [Jf(1,:)]=calcFlux(X(:,2),D,i);
        [Jf(2,:)]=calcFlux(X(:,2),D,i+1);

        Jx(i,i-1) = - Jf(1,1);
        Jx(i,i) = Jf(2,1)-Jf(1,2);
        Jx(i,i+1) = Jf(2,2);
    end
end

function [J] = calcFlux(X,D,i)
    N = size(X);
    J =[0,0];
    if and(i~=1,i~=N+1)
        x = max(X(i),X(i-1));
        imax = eq(x,X(i-1)) + eq(x,X(i))*2;
        for j=1:2
            J( j ) = D( j ) * (2 * x * eq( imax, j ) * (X( j ) - X( j - 1 )) + x ^ 2 * (eq( j, 2 ) - eq( j, 1 )));
        end
    end
end
