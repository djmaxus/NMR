import numpy as np
import time


class TimeLog():
    def __init__(self, linsol=0, proj=0, pod=0, resbld=0, jacbld=0, Nt=0):
        self.linsol = linsol
        self.proj = proj
        self.pod = pod
        self.resbld = resbld
        self.jacbld = jacbld
        self.kn = np.zeros((Nt+1))


class Solution():
    def __init__(self, param):
        self.param = param
        self.timelog = TimeLog(Nt=param.Nt)
        self.is_encoded = False
        self.is_gappy = False
        self.X = np.zeros((param.Nx, param.Nt+1))
        self.t = np.zeros((param.Nt + 1))
        self.Rx = np.zeros((param.Nx + 1, 1))
        self.Jnz = np.zeros((param.Nx+1,1))


class Diffusion():
    def __init__(self, Nx, Amp=0, Period=0, Scale=1, model=None):
        self.val = np.zeros((Nx+1, 1))
        for i in range(Nx+1):
            self.val[i] = (
                Amp * np.sin(Period * 2*np.pi * (i+1) / (Nx+1)) + 1) * Scale
        self.model = model


class Boundary():
    def __init__(self, left, right, kind=1):
        self.first = (kind == 1)
        self.val = np.array((left, right))


def setInitial(Nx, Amp=0, Period=0, Scale=1):
    x0 = np.zeros((Nx, 1))
    for i in range(Nx):
        x0[i] = Amp * np.sin(Period * 2*np.pi * (i+1) / Nx) + Scale
    return x0


def setSource(Nx, Pos=1.0, Val=0):
    q = np.zeros((Nx, 1))
    index = round(Nx*Pos)
    q[index-1] = Val
    return q


class Parameters():
    def __init__(self, Nx=None, Nt=None, D=None, bc=None, x0=None, src=None):
        self.Nx = Nx
        self.Nt = Nt
        self.D = D
        self.bc = bc
        self.x0 = x0
        self.q = src


def setRomStage(rom_stage):
    is_encoded = rom_stage >= 1
    is_gappy = rom_stage >= 2
    return is_encoded, is_gappy


def calcFluxRes(X, param, i):
    val = 0
    if param.bc.first:
        delta = np.zeros((2, 1))
        if i == 0:
            delta[0] = X[i]
            delta[1] = param.bc.val[0]
        elif i == param.Nx:
            delta[0] = param.bc.val[1]
            delta[1] = X[i - 1]
        else:
            delta = np.array([X[i], X[i - 1]])
        x = max(delta[0], delta[1])
        mult = param.D.model(x)
        val = - param.D.val[i]*param.Nx**2 * mult[0] * (delta[0] - delta[1])
    else:
        pass

    return val


def buildResidual(X, dt, param, encoder=None):
    Rt = (X[:, 1] - X[:, 0]) / dt
    Rt = np.reshape(Rt,   (param.Nx,   1))
    mask = range(param.Nx+1)
    W = np.zeros((param.Nx + 1, 1))
    for i in mask:
        W[i] = calcFluxRes(X[:, 1], param, i)
    return Rt, W


def calcFluxJac(X, param, i):
    val = np.zeros((1, 2))
    if param.bc.first:
        delta = np.zeros((2, 1))
        if i == 0:
            delta[0] = param.bc.val[1]
            delta[1] = X[i]
        elif i == param.Nx:
            delta[0] = X[i-1]
            delta[1] = param.bc.val[0]
        else:
            delta[0] = X[i-1]
            delta[1] = X[i]
        imax = 0
        if delta[1] >= delta[0]:
            imax = 1
        x = delta[imax]
        mult = param.D.model(x)
        for j in range(2):
            x_dx = 0
            if imax == j:
                x_dx = mult[1]

            d_dx = -1
            if j == 1:
                d_dx = 1

            if i == 0:
                if imax == 0:
                    x_dx = 0
                if j == 0:
                    d_dx = 0
            elif i == param.Nx:
                if imax == 1:
                    x_dx = 0
                if j == 1:
                    d_dx = 0

            val[0,j] = -param.D.val[i] * param.Nx**2*(x_dx *   (delta[1]-delta[0])+mult[0]*d_dx)
    else:
        pass
    return val


def buildFlow(X, dt, param, encoder=None):
    Jt = np.eye(param.Nx)/dt
    Jf = np.zeros((param.Nx + 1, 2))
    mask = range(param.Nx + 1)
    for i in mask:
        Jf[i, :] = calcFluxJac(X[:, 1], param, i)
    return Jt, Jf


def buildJacobian(Jf, Nx):
    Jx = np.zeros((Nx, Nx))
    for i in range(Nx):
        if i != 0:
            J1 = Jf[i, 0]
            Jx[i, i-1] = -J1
        J1 = Jf[i, 1]
        J2 = Jf[i + 1, 0]
        Jx[i, i] = J2 - J1
        if i + 1 != Nx:
            J2 = Jf[i+1, 1]
            Jx[i, i+1] = J2
    return Jx


def solve(param, rom_stage=0, encoder=None, train_mode=False, adaptive_fom=False, is_step_fom=False):
    return_code = 0
    return_message = 'OK'
    # Parameter save
    sol = Solution(param)
    # encoding status
    is_encoded, is_gappy = setRomStage(rom_stage)
    sol.is_encoded = is_encoded
    sol.is_gappy = is_gappy

    # time settings
    t = 0.0
    dt = 1./param.Nt
    dt_min = dt
    dt_max = dt

    # initial condition
    X = np.concatenate((param.x0, param.x0), axis=1)
    if is_encoded:
        Y = encoder.x.encode(X)

    # Newton iteration parameters
    kmax = 12
    kmin = 4
    crit_rel = 1e-1
    crit_abs = 1e-3

    # dataset initialization
    #     if is_encoded
    #         sol.Y = zeros(POD.X.npc(),round(1/dt_min)+1);
    #     end

    # solution process
    R0 = 1
    nstep = 0
    while t < 1.0:
        dt = min(dt, 1 - t)
        if is_encoded and (not is_step_fom):
            Y[:, 1] = Y[:, 0]
        else:
            X[:, 1] = X[:, 0]

        sol.t[nstep] = t
        sol.X[:, nstep] = X[:, 0]

        is_converged = False
        for k in range(kmax):
            if is_encoded and (not is_step_fom):
                t_pod = -time.time()
                X[:, 1] = encoder.x.decode(Y[:, 1])
                t_pod += time.time()
                sol.timelog.pod += t_pod

            # Residual
            t_res = -time.time()
            Rt, W = buildResidual(X, dt, param,  encoder)
            t_res += time.time()
            sol.timelog.resbld += t_res
            Rx = W[1:param.Nx+1] - W[0:param.Nx]
            Rs = param.q
            R = Rt + Rx + Rs
            # Jacobian
            t_jac = -time.time()
            Jt, Jf = buildFlow(X,  dt,  param,  encoder)
            t_jac += time.time()
            sol.timelog.jacbld += t_jac
            Jx = buildJacobian(Jf,  param.Nx)
            J = Jt + Jx

            if train_mode:
                sol.Rx = np.append(sol.Rx, W, axis=1)
                sol.Jnz = np.append(sol.Jnz, np.reshape(Jf[:, 0],(param.Nx+1,1)), axis=1)
                sol.Jnz = np.append(sol.Jnz,-np.reshape(Jf[:, 1],(param.Nx+1,1)),axis = 1)

            # Decoder Jacobian
            if is_encoded and (not is_step_fom):
                t_proj = -time.time()
                J = J @ encoder.x.decJac(Y[:,  1])
                t_proj += time.time()
                sol.timelog.proj += t_proj

            # convergence
            delta = np.linalg.norm(J.T@R)
            if k == 0:
                R0 = delta

            is_conv_abs = (delta * dt) <= crit_abs
            is_conv_rel = delta <= (crit_rel*R0)
            is_converged = is_conv_abs or is_conv_rel

            if is_converged:
                sol.timelog.kn[nstep-1] = k
                nstep += 1
                t += dt

                if is_encoded and (not is_step_fom):
                    Y[:,  0] = Y[:,  1]
                    X[:,  0] = encoder.x.decode(Y[:,  0])
                else:
                    X[:, 0] = X[:, 1]

                # timestep increase
                if k < kmin:
                    dt *= 2
                    dt = min(dt,   dt_max)

                if is_encoded and is_step_fom:
                    Y = encoder.x.encode(X)

                break
            elif (k + 1 == kmax):
                dt /= 4
                break

            # linear system solve
            t_lin = -time.time()
            dX = np.linalg.solve(J.T @ J,   -J.T@R)
            t_lin += time.time()
            sol.timelog.linsol += t_lin

            is_weak = np.sum(np.isnan(dX)) + np.sum(np.isinf(dX))
            if is_weak > 0:
                dX[:] = 0

            # update
            if is_encoded and (not is_step_fom):
                Y[:, 1] += dX
            else:
                X[:, 1] += np.reshape(dX,   param.Nx)

        if is_converged:
            is_step_fom = not (is_step_fom and (dt >= 1/param.Nt))
        elif (k+1) >= kmax:
            if adaptive_fom and (not is_step_fom):
                is_step_fom = True
                continue

            if dt < dt_min:
                return_message = 'Not converged'
                break
            dt /= 4

    # sol.t(nstep) = t;
    # sol.t(nstep+1:end) = [];
    # sol.X(:,nstep) = X(:,1);
    # sol.X(:,nstep+1:end) = [];
    # if is_encoded
    #    sol.Y(:,nstep) = Y(:,1);
    #    sol.Y(:,nstep+1:end) = [];
    # end
    #sol.timelog.bld = sol.timelog.resbld + sol.timelog.jacbld + sol.timelog.pod;
    #sol.timelog.lin = sol.timelog.linsol + sol.timelog.proj;
    #sol.timelog.overall = sol.timelog.lin + sol.timelog.bld;
    #sol.timelog.iters.bld = sum(sol.timelog.kn,'all');
    #sol.timelog.iters.lin = sum(sol.timelog.kn-1,'all');
    #sol.timelog.lin_avg = sol.timelog.lin/sol.timelog.iters.lin;
    #sol.timelog.bld_avg = sol.timelog.bld/sol.timelog.iters.bld;
    #sol.timelog.iters.lin_avg = sol.timelog.iters.lin / length(sol.timelog.kn);
    #sol.timelog.iters.bld_avg = sol.timelog.iters.bld / (length(sol.timelog.kn)-1);
    return sol, return_code, return_message
