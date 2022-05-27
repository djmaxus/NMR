function [out] = calcSpeedUp(sol, ref)
%CALCSPEEDUP Summary of this function goes here
%   Detailed explanation goes here
    out.sol = calcTime(sol);
    out.ref = calcTime(ref);
    
    out.rel = calcRel(out.sol,out.ref);
end

function out = calcTime(sol)
    tbl = struct2table(sol.nlin);
    
    out.build = sum(tbl.build(:));
    out.solve = sum(tbl.solve(:)); 
    out.proj  = sum(tbl.proj(:)); 
    
    out.all = out.build + out.solve + out.proj;
end

function out = calcRel(sol,ref)
out.build = ref.build / sol.build;
out.solve = ref.solve / sol.solve;
out.proj  = ref.proj  / sol.proj;
out.all   = ref.all   / sol.all;
end

