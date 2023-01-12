classdef MRST_FOM
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
        var_pod;
        res_pod;
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
            [data, res] = simulateScheduleAD(obj.state0, model,schedule);
        end
        function obj = makeSnapshots(obj, param_points)
            out = cell(length(param_points),1);
            for i = 1:length(param_points)
                out{i} = struct();
                [out{i}.res_state, out{i}.well_state] = obj.solve(param_points{i});
            end
            obj.snapshots = out;
        end
        function obj = makeVarPOD(obj,energy)
            % build sample matrix
            var_data = [];
            for i=1:obj.snapshots
                % get field var vector
                var_evol = [];
                for j = 1:length(obj.snapshots{i})
                    var_step = obj.snapshots{i}.pres;
                    var_evol = [var_evol,var_step];
                end
                % add to matrix
                var_data = [var_data,var_evol];
            end
            % make POD
            obj.var_pod = makePOD(var_data,energy);
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
