function [Jx] = buildJacobian(Jf,N)
    Jx = spalloc(N,N,3*N-2);
    for i=1:N
        if i ~=1
            J1 = Jf(i,1);
            Jx(i,i-1) = -J1;
        end
        J1 = Jf(i,2);
        J2 = Jf(i+1,1);
        Jx(i,i) = J2 - J1;
        if i ~= N
            J2 = Jf(i+1,2);
            Jx(i,i+1) = J2;
        end
    end
end
