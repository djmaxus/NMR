param_gen_h = 2^(-3);
param_gen_clx = 2^(-3);

gen_initial = @(Nu) genField(Nu, 0, 2, param_gen_h, param_gen_clx);

gen_q_pos = @(k) (-4)^(k) * (rand()-0.5)^(2*k+1)+0.5;

gen_kappa = @(Nf) genField(Nf, 2e-3, 2e-4, param_gen_h, param_gen_clx);

Nf = 2^11;
Nu = Nf + 1;

fom = FOM(Nu,1,1,gen_initial, gen_kappa,struct('pos',gen_q_pos(2),'value',1.0));

fom.plot_problem();

Nt = Nf;

tic;
[sol] = ode_solve(@fom.calc_linearization,fom.u0,Nt,true);
toc

fom.plot_solution(); 

%% POD builder

builder_pod = BuilderPOD(sol.U, fom.u0);

%% POd encoder

[builder_pod,encoder_pod] = builder_pod.set_energy_residual(1e-3).build_encoder();
