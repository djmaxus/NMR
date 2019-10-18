function [wellState, resState] = ResSimSPCA(problem)
%RESSIMSPCA Summary of this function goes here
%   Detailed explanation goes here

%% init computational grid
G = cartGrid(problem.grid.N, problem.grid.L * meter());
G = computeGeometry(G);

%% init rock
rock = makeRock(G,problem.rock.perm * milli * darcy, problem.rock.poro);

%% init fluid
fluid = initSimpleADIFluid('phases',problem.fluid.phases,...
    'mu',problem.fluid.mu * centi() *poise(),...
    'rho',problem.fluid.rho * kilogram() / meter()^3,...
    'c',problem.fluid.c / barsa(),...
    'cR',problem.rock.cR / barsa(),...
    'pRef',problem.pRef * barsa());

%% init model
model = TwoPhaseOilWaterModel(G,rock,fluid);

%% init wells
well = struct([]);
for w = 1:length(problem.well)
    obj = problem.well(w);
    pos = obj.pos;
    index = sub2ind(G.cartDims,pos(1),pos(2),pos(3));
    well = addWell(well,G,rock,index,...
        'Type',obj.type,...
        'Val', obj.val * barsa(),...
        'Radius',obj.radius * meter(),...
        'Sign',2 * obj.is_inj - 1,...
        'Comp_i',obj.comp);
end


%% set boundary condition
bc = struct([]);
bc = fluxside(bc,G,'left',0,'sat',[0,0]);
bc = fluxside(bc,G,'right',0,'sat',[0,0]);
bc = fluxside(bc,G,'back',0,'sat',[0,0]);
bc = fluxside(bc,G,'front',0,'sat',[0,0]);
bc = fluxside(bc,G,'bottom',0,'sat',[0,0]);
bc = fluxside(bc,G,'top',0,'sat',[0,0]);


%% init time
T = problem.time.T * year();
dt = rampupTimesteps(T,problem.time.step*day(),0);

%% init schedule
schedule = simpleSchedule(dt,'W',well,'bc',bc);

%% set initial condition
P0 = problem.init.P * barsa();
S0 = [1 - problem.init.S, problem.init.S];
state0 = initResSol(G,P0,S0);

%% set up output handling

%% solve problem
[wellState,resState] = simulateScheduleAD(state0,model,schedule);
%% return output
% plotCellData(G, resState{end}.s(:,1));
end
