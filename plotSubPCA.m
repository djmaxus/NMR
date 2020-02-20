U = POD.U * POD.S;
hold on;
plot3([0,U(idx(1),1)],[0,U(idx(2),1)],[0,U(idx(3),1)],'-*r');
hold on;
plot3([0,U(idx(1),2)],[0,U(idx(2),2)],[0,U(idx(3),2)],'-*r');
hold on;
plot3([0,U(idx(1),3)],[0,U(idx(2),3)],[0,U(idx(3),3)],'-*r');
