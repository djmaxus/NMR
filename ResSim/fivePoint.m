function [problem] = fivePoint(perm)
    %% variable parameters
    cell_num = length(perm);
    Nx = sqrt(cell_num);
    N = [Nx;Nx;1];
    perm = reshape(perm,[cell_num,1]);
    %% grid
    grid.N = N;
    grid.L = [410;410;3];
    problem.grid = grid;
    %% rock
    rock.perm = perm;
    rock.poro = 0.3;
    rock.cR = 4.35e-5;
    
    problem.rock = rock;
    
    %% pRef
    problem.pRef = 90;
    
    %% fluid
    fluid.phases = 'OW';
    fluid.mu = [2.2;1.0];
    fluid.rho = [800;1000];
    fluid.c = [8e-5;4e-5];
    fluid.B = [1.0;1.0];
    
    problem.fluid = fluid;
    
    %% well
    well_base.type = 'bhp';
    well_base.radius = 0.07;
    well_base.pos = [0,0,1];
    
    prod = well_base;
    prod.val = 100.0;
    prod.is_inj = false;
    prod.comp = [1;0;0];
    
    inj = well_base;
    inj.val = 350;
    inj.is_inj = true;
    inj.comp = [0;1;0];
    inj.pos = [round((grid.N(1)+1)/2),round((grid.N(2)+1)/2),1];
    
    well = struct(prod);
    
    well(1:4) = prod;
    wIndX = [1,grid.N(1)];
    wIndY = [1,grid.N(2)];
    prod_val = [100,110,120,130];
    for i = 1:2
        for j = 1:2
            well(2*(i-1) + j).pos = [wIndX(i), wIndY(j),1];
            well(2*(i-1) + j).val = prod_val(2*(i-1) + j);
        end
    end
    well(5) = inj;
    problem.well = well;
    %% time
    time.T = 10;
    time.step = 60;
    
    problem.time = time;
    
    %% init
    init.P = 100;
    init.S = 0.0;
    
    problem.init = init;
    
end

