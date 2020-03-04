nR = 0.1;
nR = round(nR * Nx);
nR = max(nR,POD.Rx.N);
POD.Rx = makeOblique(nR/Nx,POD.Rx);
disp(norm(POD.Rx.O * data_r - data_r));

nJ = 0.1;
nJ = round(nJ * size(data_j,1));
nJ = max(nJ,POD.Jnz.N);
POD.Jnz = makeOblique(nJ/Nx/3,POD.Jnz);
disp(norm(POD.Jnz.O * data_j - data_j));
