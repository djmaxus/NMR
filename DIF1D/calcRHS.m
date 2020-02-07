function [R,J] = calcRHS(X,dt,D,dx,bc)
    N = size(X,1);
    R = zeros(N,1);
    J = spalloc(N,N,3*N-2);
    % time derivative
    for i=1:N
        R(i) = (X(i,2)-X(i,1))/dt;
        J(i,i) = 1/dt;
    end
    % divergence
    flux(1) = bc(1);
    for i=1:N
        if i==N
            flux(2)=bc(2);
            Jf(2,:) = 0;
        else
            [flux(2),Jf(2,:)]=calcFlux(X(:,2),D,dx,i+1);
        end

        R(i) = R(i) + flux(2) - flux(1);

        if i~=1
            J(i,i-1) = - Jf(1,1);
        end
        J(i,i) = J(i,i) + Jf(2,1)-Jf(1,2);
        if i~=N
            J(i,i+1) = Jf(2,2);
        end

        flux(1)=flux(2);
        Jf(1,:) = Jf(2,:);
    end

end

function [R,J] = calcFlux(X,D,dx,i)
    x = max(X(i),X(i-1));
    imax = eq(x,X(i-1)) + eq(x,X(i))*2;
    R = D( i )/dx(i) * x^2 * ( X( i ) - X( i - 1 ) );
    J =[0,0];
    for j=1:2
        J(j) = D(j)/dx(j)*(2*x*eq(imax,1)*( X( j ) - X( j - 1 ) )+x^2 * (eq(j,2)-eq(j,1));
    end
end
