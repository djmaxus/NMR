function [q] = genSource(Nx,pos,val)
    %GENSOURCE Summary of this function goes here
    %   Detailed explanation goes here
    q.val = spalloc(Nx,1,1);
    q.pos = round(pos*Nx);
    q.val(q.pos) = val;
end

