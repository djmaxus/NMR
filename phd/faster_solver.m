param_gen_h = 2^(-3);
param_gen_clx = 2^(-3);

gen_initial = @(Nu) genField(Nu, 0, 2, param_gen_h, param_gen_clx);

gen_q_pos = @(k) (-4)^(k) * (rand()-0.5)^(2*k+1)+0.5;

gen_kappa = @(Nf) genField(Nf, 2e-3, 2e-4, param_gen_h, param_gen_clx);

Nf = 2^10;
Nu = Nf + 1;

calc_g_all = @(u) [1 + 0.5 * u.^2, u];

calc_g = @(u, elem) (1 + 0.5 * u(elem)^2);
calc_dg = @(u) u;

calc_g_full = @(u) eval_elem_wise(calc_g,u,1:Nu);

fom = FOM(Nu,1,1,gen_initial, gen_kappa,struct('pos',gen_q_pos(2),'value',1.0), calc_g_all);

fom.plot_problem();

Nt = Nf;

%%
tic;
[sol] = ode_solve(@fom.calc_linearization,fom.u0,Nt,true);
toc

fom.plot_solution();

%% POD builder
tic;
builder_pod = BuilderPOD(sol.U, fom.u0);
toc

%% POD encoder

encoder_pod = builder_pod.set_energy_residual(1e-3).build_encoder();

%% ROM-I

tic;
[sol_encoded] = ode_solve_encoded(@fom.calc_linearization,fom.u0,Nt,true,encoder_pod);
toc

%% Nonlinear data & Gappy builder

G = zeros(size(sol_encoded.data{1},1),numel(sol_encoded.data));
dG = G;

for i = 1:numel(sol_encoded.data)
    G(:,i) = sol_encoded.data{i}(:,1);
    dG(:,i) = sol_encoded.data{i}(:,2);
end

tic;
builder_gappy = BuilderPOD(G, 0);
builder_gappy_dg = BuilderPOD(dG, 0);
toc

%% Setup Gappy POD

for builder = [builder_gappy,builder_gappy_dg]
    builder...
        .set_energy_residual(1e-3)...
        .set_lower_bound(encoder_pod.dim_latent()*2)...
        .set_gappy_mask(2)...
        .set_mask_extension_normalized([0,1])...
        .set_mask_extension(find(fom.q));
end

%% Setup Gappy Decoder

encoder_gappy = builder_gappy.build_gappy();

calc_g_full_gappy = @(u) eval_elem_wise(calc_g,u,encoder_gappy.gappy_mask);

calc_g_gappy = @(u) encoder_gappy.eval_gappy(calc_g_full_gappy,u);

%%

u = sol_encoded.U(:,end);

t = 0;
t_gappy = 0;
wb = waitbar(0);
n_repeat = 1;
for i = 1:n_repeat
    tic;
    g = calc_g_full(u);
    t = toc() + t;
    
    tic;
    g_gappy = calc_g_gappy(u);
    t_gappy = toc() + t_gappy;
    
    wb = waitbar(i/n_repeat);
end
close(wb);

disp(t_gappy/t);

%% Gappy for derivative

