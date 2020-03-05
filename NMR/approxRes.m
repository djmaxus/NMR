nR = 0.5;
nR = round(POD.Rx.N * (1 + nR));
nR = nR/size(POD.Rx.U,1);
POD.Rx = makeOblique(nR,POD.Rx);
disp(norm(applyOblique(data_r,POD.Rx) - data_r));

nJ = 0.5;
nJ = round(POD.Jnz.N * (1 + nJ));
nJ = nJ/size(POD.Jnz.U,1);
POD.Jnz = makeOblique(nJ,POD.Jnz);
disp(norm(applyOblique(data_j,POD.Jnz) - data_j));
