function [Jt,Jx,Jnz] = buildJacobian(X,dt,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);
    % time derivative
    Jt = speye(N)./dt;
    % flux 
    
    Jnz = zeros(3*N-2,1);
    mask = 1:3*N-2;
    if nmr
        mask = POD.Jnz.idx_mask';
    end
    for k = mask
        if k <= N
            i = k;
            j = k;
        elseif k <= 2*N-1
            i = k - N;
            j = i+1;
        else
            j = k - 2*N+1;
            i = j+1;
        end
        if i == j
            J1 = calcFlux(X(:,2),D,bc,i,2);
            J2 = calcFlux(X(:,2),D,bc,i+1,1);
            Jnz(k) = J2 - J1;
        elseif i < j
            J2 = calcFlux(X(:,2),D,bc,i+1,2);
            Jnz(k) = J2;
        else
            J1 = calcFlux(X(:,2),D,bc,i,1);
            Jnz(k) = -J1;
        end
    end
    if nmr
        % gappy;
        Jnz = applyOblique(Jnz,POD.Jnz);
    end
    
    Jc = Jnz(1:N);
    Ju = Jnz(N+1:2*N-1);
    Jd = Jnz(2*N:3*N-2);
    
    Jx = spalloc(N,N,3*N-2);
    for i=1:N
        if i ~=1
            Jx(i,i-1) = Jd(i-1);
        end
        Jx(i,i) = Jc(i);
        if i ~= N
            Jx(i,i+1) = Ju(i);
        end
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
        for j = rj
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
            for j=rj
                J( j ) = -D( i ) * (2 * x * eq( imax, j ) * (X( i ) - X( i - 1 )) + (x^2+1) * (eq( j, 2 ) - eq( j, 1 )));
            end
        end
    end
    J = J(rj);
end
