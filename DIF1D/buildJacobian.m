function [Jt,Jx] = buildJacobian(X,dt,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);
    % time derivative
    Jt = ones(N,1)./dt;

    % flux
    Jx = spalloc(N,N,3*N-2);
    %TODO: make 3 diags
    mask = 1:N;
    if nmr
        %mask = POD.Rx.idx_mask';
    end
    Jw = zeros(N+1,2);
    for i = mask
        Jw(i,:) = calcFlux(X(:,2),D,bc,i);
        Jw(i+1,:) = calcFlux(X(:,2),D,bc,i+1);
    end
    for i=1:N
        if i ~=1
            Jx(i,i-1) = -Jw(i,1);
        end
        Jx(i,i) = Jw(i+1,1) - Jw(i,2);
        if i ~= N
            Jx(i,i+1) = Jw(i+1,2);
        end
    end
    if nmr
        Jx = applyOblique(Jx',POD.Rx)';
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
