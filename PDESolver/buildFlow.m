function [Jt,Jf] = buildFlow(X,dt,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);
    % time derivative
    Jt = speye(N)./dt;
    % flux
    Jf = zeros(N+1,2);
    mask = 1:N+1;
    if nmr
        mask = POD.idx_mask;
    end
    for i = mask
        Jf(i,:) = calcFlux(X(:,2),D,bc,i,1:2);
    end
    if nmr
        % gappy;
        Jf = applyOblique(Jf.*[1,-1],POD).*[1,-1];
    end
end

function [J] = calcFlux(X,D,bc,i,rj)
    N = length(X);
    J = zeros(2,1);
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
        mult = D.model(x);%[x^2+1,2*x];
        for j = rj
            x_dx = mult(2)*eq(imax,j);
            d_dx = eq( j, 2 ) - eq( j, 1 );
            if i == 1
                x_dx = x_dx * ne(imax,1);
                d_dx = d_dx * ne(j,1);
            elseif i == N+1
                x_dx = x_dx * ne(imax,2);
                d_dx = d_dx * ne(j,2);
            end
            J( j ) = -D.val( i ) * (x_dx * (delta( 2 ) - delta(1)) + mult(1) * d_dx);
        end
    else
        if and(i~=1,i~=N+1)
            x = max(X(i),X(i-1));
            imax = eq(x,X(i-1)) + eq(x,X(i))*2;
            for j=rj
                J( j ) = -D.val( i ) * (2 * x * eq( imax, j ) * (X( i ) - X( i - 1 )) + (x^2+1) * (eq( j, 2 ) - eq( j, 1 )));
            end
        end
    end
    J = J(rj);
end
