function [Jt,Jx,Jw] = buildJacobian(X,dt,D,bc,POD)
    N = size(X,1);
    nmr = ~isempty(POD);
    % time derivative
    Jt = ones(N,1)./dt;
    % flux 
    Jw = zeros(N+1,2);
    mask = 1:N+1;
    if nmr
        mask = POD.Jw.idx_mask';
    end
    for k = mask
        % set i,j
        rj = 1:2;
        i = k;
        if k > N + 1
            i = k - N - 1;
            rj = 2;
        end
        Jw(i,rj) = calcFlux(X(:,2),D,bc,i,rj);
    end
    if nmr
        % gappy;
        Jw = [Jw(:,1);Jw(:,2)];
        Jw = applyOblique(Jw,POD.Jw);
        Jw = [Jw(1:N+1,1),Jw(N+2:end,1)];
    end
    
    Jc = zeros(N,1);
    Ju = zeros(N-1,1);
    Jd = zeros(N-1,1);
    for i=1:N
        if i ~=1
            Jd(i-1) = -Jw(i,1);
        end
        Jc(i) = Jw(i+1,1) - Jw(i,2);
        if i ~= N
            Ju(i) = Jw(i+1,2);
        end
    end
    Jw = [Jw(:,1);Jw(:,2)];
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
    J = zeros(size(rj));
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
