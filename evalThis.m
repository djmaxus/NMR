function out = evalThis(input)
%  EVALTHIS Brief summary of this function.
% 
% Detailed explanation of this function.
    data = plot(axes(figure),input,input.^2);
    
    save('buf.mat');
    out = load('buf.mat');
    delete("buf.mat");
end