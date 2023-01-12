function [is_type_correct] = check_class(instance,reference)
%CHECK_CLASS Summary of this function goes here
%   Detailed explanation goes here
sc = superclasses(instance);
sc = {sc,class(instance)};

is_reference_found = false;
for cl = sc
    if strcmp(cl{:},reference)
        is_reference_found = true;
        break;
    end
end
is_type_correct = is_reference_found;
if ~is_type_correct
    error('Incorrect input type');
end
end

