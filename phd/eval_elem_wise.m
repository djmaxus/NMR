function [out] = eval_elem_wise(func,arg,mask)
arguments
    func
    arg
    mask (1,:) {mustBeInteger, mustBePositive}
end

out = zeros(length(mask),1);

for index = 1:length(mask)
    out(index) = func(arg, mask(index));
end

end
