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
    for i=1:length(mask)
        j = mask(i);
        Jf = [0,0;0,0];
        [Jf(1,:)]=calcFlux(X(:,2),D,j);
        [Jf(2,:)]=calcFlux(X(:,2),D,j+1);
        if j ~=1
            Jx(j,j-1) = - Jf(1,1);
        end
        
        Jx(j,j) = Jf(2,1)-Jf(1,2);
        if j ~= N
            Jx(j,j+1) = Jf(2,2);
        end
    end
end

function [J] = calcFlux(X,D,i)
    N = size(X);
    J =[0,0];
    if and(i~=1,i~=N+1)
        x = max(X(i),X(i-1));
        imax = eq(x,X(i-1)) + eq(x,X(i))*2;
        for j=1:2
            J( j ) = -D( i ) * (2 * x * eq( imax, j ) * (X( i ) - X( i - 1 )) + x ^ 2 * (eq( j, 2 ) - eq( j, 1 )));
        end
    end
end
