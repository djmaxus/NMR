%% variable snapshots
trainVar;
POD.X = makePOD(data_x,true,false,0.99);
POD.X.data = cell(size(vals,1),1);
for i = 1:length(POD.X.data)
    POD.X.data{i} = data{i};
end
%% residual snapshots
trainRes;

%% hyper reduction
truncRes;
approxRes;
