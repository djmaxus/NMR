classdef FOM
    %COSTOBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        well_sol;
        res_sol;
        problem;
        G;
        fluid;
        bc;
        dt;
        state0;
        snapshots;
    end
    
    methods
        function obj = FOM(X)
            %COSTOBJECT Construct an instance of this class
            %   Detailed explanation goes here
            problem = fivePoint(X);
            obj.problem = problem;
            %% init computational grid
            obj.G = cartGrid(problem.grid.N, problem.grid.L * meter());
            obj.G = computeGeometry(obj.G);
            
            %% init fluid
            obj.fluid = initSimpleADIFluid('phases',problem.fluid.phases,...
                'mu',problem.fluid.mu * centi() *poise(),...
                'rho',problem.fluid.rho * kilogram() / meter()^3,...
                'c',problem.fluid.c / barsa(),...
                'cR',problem.rock.cR / barsa(),...
                'pRef',problem.pRef * barsa()); 
            
            %% set boundary condition
            bc = struct([]);
            bc = fluxside(bc,obj.G,'left',0,'sat',[0,0]);
            bc = fluxside(bc,obj.G,'right',0,'sat',[0,0]);
            bc = fluxside(bc,obj.G,'back',0,'sat',[0,0]);
            bc = fluxside(bc,obj.G,'front',0,'sat',[0,0]);
            bc = fluxside(bc,obj.G,'bottom',0,'sat',[0,0]);
            bc = fluxside(bc,obj.G,'top',0,'sat',[0,0]);
            obj.bc = bc;
            
            %% init time
            T = problem.time.T * year();
            obj.dt = rampupTimesteps(T,problem.time.step*day(),0);
            
            %% set initial condition
            P0 = problem.init.P * barsa();
            S0 = [1 - problem.init.S, problem.init.S];
            obj.state0 = initResSol(obj.G,P0,S0);
            
            %% get hist data
            [obj.well_sol, obj.res_sol] = obj.solve(X);
        end
        function [data, res] = solve(obj,X)
            %% init rock
            rock = makeRock(obj.G,X * milli * darcy, obj.problem.rock.poro);
            %% init model
            model = TwoPhaseOilWaterModel(obj.G,rock,obj.fluid);
            %% init wells
            well = struct([]);
            for w = 1:length(obj.problem.well)
                item = obj.problem.well(w);
                pos = item.pos;
                index = sub2ind(obj.G.cartDims,pos(1),pos(2),pos(3));
                well = addWell(well,obj.G,rock,index,...
                    'Type',item.type,...
                    'Val', item.val * barsa(),...
                    'Radius',item.radius * meter(),...
                    'Sign',2 * item.is_inj - 1,...
                    'Comp_i',item.comp);
            end
            %% init schedule
            schedule = simpleSchedule(obj.dt,'W',well,'bc',obj.bc);
            disp('Simulating...');
            disp('');
            [data, res] = simulateScheduleAD(obj.state0, model,schedule);
        end
        function obj = makeSnapshots(obj)
            obj.snapshots = [];
            min_val = [obj.res_sol{1}.pressure(1);obj.res_sol{1}.s(1,2)];
            max_val = min_val;
            Nt =length(obj.res_sol);
            Nc = length(stamp.pressure);
            for t = 1:length(obj.res_sol)
                stamp = obj.res_sol{t};
                snapshot = [];
                for i = 1:length(stamp.pressure)
                    val = [stamp.pressure(i);stamp.s(i,2)];
                    max_val = max(max_val,val);
                    min_val = min(min_val,val);
                    snapshot = [snapshot;val];
                end
                obj.snapshots = [obj.snapshots,snapshot];
            end
        end
        function [U,S,V,Nr] = makePOD(obj,w)
            [U,S,V] = svd(obj.snapshots);
            energy = diag(S).^2;
            total_energy = sum(energy);
            fraction = 0;
            for i = 1:length(energy)
                fraction = fraction + energy(i)/total_energy;
                if fraction >= w
                    Nr = i;
                    break;
                end
            end
        end
        function plot_solution(obj,res)
            plotCellData(obj.G, res{end}.s(:,1));
        end
        function [data,cost] = calc(obj,X)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [data, ~] = obj.solve(X);
            cost = obj.getCost(data);
        end
        
        function [cost] = getCost(obj,data)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            cost = 0;
            Tn = length(data);
            for t = Tn
                step = data(t);
                step_hist = obj.hist_data(t);
                Wn = length(step{1});
                for w = 1:Wn
                    well = step{1}(w);
                    well_hist = step_hist{1}(w);
                    
                    flux = well.flux;
                    bhp = well.bhp;
                    flux_hist = well_hist.flux;
                    bhp_hist = well_hist.bhp;
                    
                    sol = [flux,bhp];
                    sol_hist = [flux_hist,bhp_hist];
                    Sn = length(sol);
                    denom = Tn * Wn * Sn; 
                    denom = denom + eq(denom,0);
                    weight = 1./(sol_hist.^2 .* ne(sol_hist,0) + eq(sol_hist,0));
                    delta = (sol - sol_hist).^2 .* weight ./denom;
                    cost = cost + sum(delta);
                end
            end
            disp('');
            disp(cost);
            disp('');
        end
    end
end

