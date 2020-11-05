import numpy as np
import solver

Nt = 10
Nx = 10
param = solver.Parameters(D=solver.Diffusion(val=np.ones((Nx, 1))))
param.x0 = np.zeros((Nx, 1))
sol, code, message = solver.solve(Nx, Nt, param)
print(message)
