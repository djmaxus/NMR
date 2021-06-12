function [R] = calcFlowRes(X,D,bc,i)
    N = length(X);
    if bc.kind == 1
        if i==1
            delta = [X( i ),bc.state(1,2)];
        elseif i==N+1
            delta = [bc.state(2,2), X( i - 1 )];
        else
            delta = [X( i ),X( i - 1)];
        end
        x = max(delta(2),delta(1));
        mult = D.model(x);%(x^2+1);
        R = -D.val( i ) * mult(1) * ( delta(1) - delta(2) );
    else
        if i==1
            R = bc.val(1);
        elseif i==N+1
            R = bc.val(2);
        else
            x = max(X(i),X(i-1));
            R = -D.val( i ) * (x^2+1) * ( X( i ) - X( i - 1 ) );
        end
    end
    R = R * (N^2);
end

