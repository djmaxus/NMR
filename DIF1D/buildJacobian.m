function [Jt,Jx] = buildJacobian(X,dt,D,bc,mask,POD)
    N = size(X,1);
    nmr = length(mask) < N;
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
        [Jf(1,:)]=calcFlux(X(:,2),D,bc,j);
        [Jf(2,:)]=calcFlux(X(:,2),D,bc,j+1);
        if j ~=1
            Jx(j,j-1) = - Jf(1,1);
        end
        
        Jx(j,j) = Jf(2,1)-Jf(1,2);
        if j ~= N
            Jx(j,j+1) = Jf(2,2);
        end
    end
    if nmr
        Jx = POD.Rx.Ut*POD.Rx.O*Jx;
    end
end

function [J] = calcFlux(X,D,bc,i)
    N = length(X);
    J =[0,0];
    if bc.first
        if i==1
            delta = [bc.val(1),X(i)];
        elseif i==N+1
            delta = [X( i - 1 ),bc.val(2)];
        else
            delta = [X( i-1 ),X( i )];
        end
        imax = 1+ ge(delta(2)-delta(1),0);
        x = delta(imax);
        for j=1:2
            x_dx = 2*x*eq(imax,j);
            d_dx = eq( j, 2 ) - eq( j, 1 );
            if i == 1
                x_dx = x_dx * ne(imax,1);
                d_dx = d_dx * ne(j,1);
            elseif i == N+1
                x_dx = x_dx * ne(imax,N+1);
                d_dx = d_dx * ne(j,2);
            end
            J( j ) = -D( i ) * (x_dx * (delta( 2 ) - delta(1)) + (x^2+1) * d_dx);
        end
    else
        if and(i~=1,i~=N+1)
            x = max(X(i),X(i-1));
            imax = eq(x,X(i-1)) + eq(x,X(i))*2;
            for j=1:2
                J( j ) = -D( i ) * (2 * x * eq( imax, j ) * (X( i ) - X( i - 1 )) + (x^2+1) * (eq( j, 2 ) - eq( j, 1 )));
            end
        end
    end
end
