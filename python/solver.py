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
        self.t = np.zeros((param.Nt+1))


class Diffusion():
    def __init__(self, val=None, model=None):
        self.val = val
        self.model = model


class Parameters():
    def __init__(self, Nx=None, Nt=None, D=None, bc=None, x0=None, q=None):
        self.Nx = Nx
        self.Nt = Nt
        self.D = D
        self.bc = bc
        self.x0 = x0
        self.q = q


def setRomStage(rom_stage):
    is_encoded = rom_stage >= 1
    is_gappy = rom_stage >= 2
    return is_encoded, is_gappy


def buildResidual(X, dt, param, encoder=None):
    Rt = np.zeros((param.Nx,   1))
    W = np.zeros((param.Nx+1,   1))
    return Rt, W


def buildFlow(X, dt, param, encoder=None):
    Jt = np.zeros((param.Nx,   param.Nx))
    Jf = np.zeros((param.Nx,   2))
    return Jt, Jf


def buildJacobian(Jf, Nx):
    Jx = np.zeros((Nx,  Nx))
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
    x0 = param.x0
    X = np.concatenate((x0, x0), axis=1)
    if is_encoded:
        Y = encoder.x.encode(X)

    # paramter distribution
    D = param.D
    D.val = (param.Nx ** 2) * D.val

    # Newton iteration parameters
    kmax = 12
    kmin = 4
    crit_rel = 1e-4
    crit_abs = 1e-4

    # dataset initialization
    #     if is_encoded
    #         sol.Y = zeros(POD.X.npc(),round(1/dt_min)+1);
    #     end
    #     if train_mode
    #         sol.Rx = [];
    #         sol.Jnz = [];
    #     end

    # solution process
    R0 = 1
    nstep = 1
    while t < 1.0:
        dt = min(dt, 1 - t)
        t += dt
        if is_encoded and (not is_step_fom):
            Y[:, 1] = Y[:, 0]
        else:
            X[:, 1] = X[:, 0]

        sol.t[nstep-1] = t
        sol.X[:, nstep - 1] = X[:, 0]
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

            #             if train_mode
            #                 % collect residual snapshots
            #                 sol.Rx = [sol.Rx,W];
            #                 sol.Jnz = [sol.Jnz,Jf.*[1,-1]];
            #             end

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
            #             if and(is_encoded,~is_step_fom)
            #                 Y(:,2) = Y(:,2) + dX;
            #             else
            #                 X(:,2) = X(:,2) + dX;
            #             end
            #         end

            #         if is_converged
            #             waitbar(t/1,wb,[wb_msg,msg_add(is_step_fom)]);
            #             if and(is_step_fom,dt >= 1/Nt)
            #                 is_step_fom = false;
            #             end
            #         end

            #         if and(~is_converged,k >= kmax)
            #             if adaptive_fom
            #                 if (~is_step_fom)
            #                     is_step_fom = true;
            #                     continue;
            #                 end
            #             end
            #             if (dt < dt_min)
            #                 warning(['Not converged',newline(),...
            #                     'Relative residual: ', num2str(delta/R0)]);
            #                 break;
            #             end
            #             dt = dt/4;
            #         end

            #     end

            #     sol.t(nstep) = t;
            #     sol.t(nstep+1:end) = [];

            #     sol.X(:,nstep) = X(:,1);
            #     sol.X(:,nstep+1:end) = [];
            #     if is_encoded
            #         sol.Y(:,nstep) = Y(:,1);
            #         sol.Y(:,nstep+1:end) = [];
            #     end

            #     sol.timelog.bld = sol.timelog.resbld + sol.timelog.jacbld + sol.timelog.pod;
            #     sol.timelog.lin = sol.timelog.linsol + sol.timelog.proj;
            #     sol.timelog.overall = sol.timelog.lin + sol.timelog.bld;
            #     sol.timelog.iters.bld = sum(sol.timelog.kn,'all');
            #     sol.timelog.iters.lin = sum(sol.timelog.kn-1,'all');
            #     sol.timelog.lin_avg = sol.timelog.lin/sol.timelog.iters.lin;
            #     sol.timelog.bld_avg = sol.timelog.bld/sol.timelog.iters.bld;
            #     sol.timelog.iters.lin_avg = sol.timelog.iters.lin / length(sol.timelog.kn);
            #     sol.timelog.iters.bld_avg = sol.timelog.iters.bld / (length(sol.timelog.kn)-1);
            #     delete(wb);
            # end
    return sol, return_code, return_message
