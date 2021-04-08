function vals = setParSample(pMax)
    vals = [];
    N = 16;
    %% limit values

    out = zeros(1,2); % zero point
    vals = [vals,out];
    N = N - 1;

    out = ones(2,2); % same limit
    out(2,:) = -1;
    vals = [vals;out];
    N = N - 2;

    out = ones(2,2); % opposite limit
    out(1,1) = -1;
    out(2,2) = -1;
    vals = [vals;out];
    N = N - 2;

    out = eye(2,2); % combined with zero
    vals = [vals;out;-out];
    N = N - 4;
    %
    %% scale
    vals = vals * pMax;

    %% inner samples
    out = lhsdesign(max(N,4),2,'smooth','off');
    out = 2 * pMax * out - pMax;
    vals = [vals;out];
end
%% plot
% f = figure;
% f.WindowStyle = 'docked';
