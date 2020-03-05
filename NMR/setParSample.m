%% limit values
% zero point
vals = zeros(1,2);

% same limit
out = ones(2,2);
out(2,:) = -1;
vals = [vals;out];

% opposite limit
out = ones(2,2);
out(1,1) = -1;
out(2,2) = -1;
vals = [vals;out];

% combined with zero
out = eye(2,2);
vals = [vals;out;-out];

% scale
vals = vals * pMax;

%% inner samples
out = lhsdesign(4,2,'criterion','correlation');
out = 2 * pMax * out - pMax;
vals = [vals;out];
