import numpy as np
import solver


def model(u):
    return np.array([u**2+1, 2*max(u, 0)])
Nt = 10
Nx = 10
x0 = solver.setInitial(Nx, 0.5, 2, 1)
bc = solver.Boundary(0, 0)
D = solver.Diffusion(Nx, 0.9, 1, 5e-2, model=model)
q = solver.setSource(Nx, 0.6, 1e+3)
param = solver.Parameters(
    Nx=Nx,
    Nt=Nt,
    D=D,
    bc=bc,
    src=q,
    x0=x0)
sol, code, message = solver.solve(param, train_mode=True)
X = sol.X
T = sol.t
W = sol.Rx
print(message)
