import numpy as np
import time


class TimeLog():
    def __init__(self, linsol=0, proj=0, pod=0, resbld=0, jacbld=0, Nt=0):
        self.linsol = linsol
        self.proj = proj
        self.pod = pod
        self.resbld = resbld
        self.jacbld = jacbld
        self.kn = np.zeros((1, Nt))


class Solution():
    def __init__(self, Nt, Nx, param):
        self.Nt = Nt
        self.Nx = Nx
        self.param = param
        self.timelog = TimeLog(Nt=Nt)
        self.is_encoded = False
        self.is_gappy = False
        self.X = np.zeros((Nx, Nt+1))
        self.t = np.zeros((1, Nt+1))


class Diffusion():
    def __init__(self, val=None, model=None):
        self.val = val
        self.model = model


class Parameters():
    def __init__(self, D=Diffusion(), bc=None, x0=None):
        self.D = D
        self.bc = bc
        self.x0 = x0


def setRomStage(rom_stage):
    is_encoded = rom_stage >= 1
    is_gappy = rom_stage >= 2
    return is_encoded, is_gappy


def buildResidual(X, dt, D, bc):
    Rt = X[:, 1]
    W = 0
    return Rt, W


def solve(Nx, Nt, param, rom_stage=0, encoder=None, train_mode=False, adaptive_fom=False, is_step_fom=False):
    return_code = 0
    return_message = 'OK'
    # Parameter save
    sol = Solution(Nt, Nx, param)
    # encoding status
    is_encoded, is_gappy = setRomStage(rom_stage)
    sol.is_encoded = is_encoded
    sol.is_gappy = is_gappy

    # time settings
    t = 0.0
    dt = 1./Nt
    dt_min = dt
    dt_max = dt

    # initial condition
    x0 = param.x0
    X = np.concatenate((x0, x0), axis=1)
    if is_encoded:
        Y = encoder.x.encode(X)

    # paramter distribution
    D = param.D
    D.val = (Nx ** 2) * D.val

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

    #     % solution process
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
            Rt, W = buildResidual(X, dt, D, param.bc)
            t_res += time.time()
            sol.timelog.resbld += t_res

            #             t_res_build = tic;
            #             [Rt,W] = buildResidual(X,dt,D,param.bc,PODr);
            #             t_res_build = toc(t_res_build);
            #             sol.timelog.resbld = sol.timelog.resbld + t_res_build;
            #             Rx = W(2:Nx+1)-W(1:Nx);
            #             R = Rt + Rx;
            #             % source term
            #             Rs = param.q.val;
            #             R = R + Rs;

            #             % Jacobian
            #             t_jac_build = tic;
            #             [Jt,Jf] = buildFlow(X,dt,D,param.bc,PODj);
            #             t_jac_build = toc(t_jac_build);
            #             sol.timelog.jacbld = sol.timelog.jacbld + t_jac_build;
            #             [Jx] = buildJacobian(Jf,Nx);
            #             J = Jx + Jt;

            #             if train_mode
            #                 % collect residual snapshots
            #                 sol.Rx = [sol.Rx,W];
            #                 sol.Jnz = [sol.Jnz,Jf.*[1,-1]];
            #             end

            #             % Decoder Jacobian
            #             if and(is_encoded,~is_step_fom)
            #                 t_projection = tic;
            #                 J = J * POD.X.JD(Y(:,2));
            #                 t_projection = toc(t_projection);
            #                 sol.timelog.proj = sol.timelog.proj + t_projection;
            #             end

            #             % convergence
            #             delta = norm(J'*R);

            #             if k == 1
            #                 R0 = delta;
            #             end
            #             is_converged_abs = delta*dt <= crit_abs;
            #             is_converged_rel = delta <= crit_rel * R0;
            #             if is_converged_abs
            #                 is_converged = is_converged_abs;
            #             else
            #                 is_converged = is_converged_rel;
            #             end

            #             if is_converged
            #                 sol.timelog.kn(nstep) = k;
            #                 nstep = nstep + 1;
            #                 t = t+dt;
            #                 if and(is_encoded,~is_step_fom)
            #                     Y(:,1) = Y(:,2);
            #                     X(:,1) = POD.X.decode(Y(:,1));
            #                 else
            #                     X(:,1) = X(:,2);
            #                 end
            #                 % timestep increase
            #                 if k <= kmin
            #                     dt =  dt * 2;
            #                     dt = min(dt,dt_max);
            #                 end
            #                 if and(is_encoded,is_step_fom)
            #                     Y = POD.X.encode(X);
            #                 end

            #                 break;
            #             elseif k == kmax
            #                 %                 dt = dt / 4;
            #                 break;
            #             end

            #             R = J'*R;
            #             J = J' * J;
            #             % linsolve
            #             t_lin_solve = tic;
            #             dX = (J) \ (-R);
            #             t_lin_solve = toc(t_lin_solve);
            #             sol.timelog.linsol = sol.timelog.linsol + t_lin_solve;

            #             is_weak = sum(isnan(dX)) + sum(isinf(dX));
            #             if is_weak > 0
            #                 dX(:) = 0;
            #             end

            #             % update
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

    # function [is_encoded,is_gappy] = setRomStage(rom_stage)
    #     is_encoded = false;
    #     is_gappy = false;

    #     if rom_stage >= 1
    #         is_encoded = true;
    #     end
    #     if rom_stage >= 2
    #         is_gappy = true;
    #     end

    # end

    # function [out] = msg_add(is_fom)
    #     if is_fom
    #         out = 'FOM';
    #     else
    #         out = 'DEF';
    #     end
    # end
