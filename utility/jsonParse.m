function [input] = jsonParse(testname)
    %JSONPARSE Summary of this function goes here
    %   Detailed explanation goes here
    fname = [testname,'.json'];
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    input = jsondecode(str);
end

