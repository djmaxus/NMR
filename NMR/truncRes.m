POD.Rx = makePOD(data_r,false,true,0.99);
POD.Rx.N = max(POD.Rx.N,POD.X.N);

POD.Jnz = makePOD(data_j,false,false,0.99);
POD.Jnz.N = max(POD.Jnz.N,POD.X.N);
