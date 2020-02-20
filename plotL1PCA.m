[Q,~] = l1pca_BF(data,POD.N,5,5,5,true);
U = Q*POD.S(1:POD.N,1:POD.N);
hold on;
plot3([0,U(idx(1),1)],[0,U(idx(2),1)],[0,U(idx(3),1)],'-+k');
hold on;
plot3([0,U(idx(1),2)],[0,U(idx(2),2)],[0,U(idx(3),2)],'-+k');
hold on;
plot3([0,U(idx(1),3)],[0,U(idx(2),3)],[0,U(idx(3),3)],'-+k');
