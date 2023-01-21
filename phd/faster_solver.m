param_gen_h = 2^(-3);
param_gen_clx = 2^(-3);

gen_initial = @(Nu) genField(Nu, 0, 2, param_gen_h, param_gen_clx);

gen_q_pos = @(k) (-4)^(k) * (rand()-0.5)^(2*k+1)+0.5;

gen_kappa = @(Nf) genField(Nf, 2e-3, 2e-4, param_gen_h, param_gen_clx);

Nf = 2^11;
Nu = Nf + 1;

calc_g = @(u, mask) {1 + 0.5 * u(mask).^2, u(mask)};

fom = FOM(Nu,1,1,gen_initial, gen_kappa,struct('pos',gen_q_pos(2),'value',1.0), @(u) calc_g(u,1:Nu));

fom.plot_problem();

Nt = round(Nf);

%%
encoder_ident = EncoderLinear(1,0);

tic;
[sol] = ode_solve_encoded(@fom.calc_linearization,fom.u0,Nt,true,encoder_ident);
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

G = zeros(size(sol_encoded.data{1}{1},1),numel(sol_encoded.data));
dG = zeros(size(sol_encoded.data{2}{1},1),numel(sol_encoded.data));

for i = 1:numel(sol_encoded.data)
    G(:,i) = sol_encoded.data{i}{1};
    dG(:,i) = sol_encoded.data{i}{2};
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

%% merge masks

common_gappy_mask = [builder_gappy.gappy_mask, builder_gappy_dg.gappy_mask];
common_gappy_mask = unique(common_gappy_mask);

for builder = [builder_gappy,builder_gappy_dg]
    builder.set_mask_extension(common_gappy_mask);
end

%% Setup Gappy encoder

encoder_gappy = builder_gappy.build_gappy();
encoder_gappy_dg = builder_gappy_dg.build_gappy();

encoder_gappy_compound = EncoderCombined({encoder_gappy, encoder_gappy_dg});

%% gappy nonlinear function

gappy_mask_compound = encoder_gappy.gappy_mask;
calc_g_gappy = @(u) calc_g(u,gappy_mask_compound);
calc_g_approx = @(u) encoder_gappy_compound.decode(calc_g_gappy(u));

rom_ii = fom;
rom_ii.calc_g = calc_g_approx;

%% ROM-II
tic;
[sol_gappy] = ode_solve_encoded(@rom_ii.calc_linearization,rom_ii.u0,Nt,true,encoder_pod);
toc

