classdef TyreModel

    properties 

        RunData % for Run Data parameters
        Sweeps % for sweep data
        splinedata
        spline_fits
        evaluation_points
        p 
        magic_formula

        P
        IA
        SA
        V
        
        FY
        FZ
        MZ
        StoredPlots

    end

    methods(Static)
        function cleanMask = filterMask(mask, minLen)
            mask = mask(:);
            d = diff([0; mask; 0]);
            s = find(d==1);
            e = find(d==-1)-1;
            cleanMask = false(size(mask));
            for i = 1:numel(s)
                if (e(i)-s(i)+1) >= minLen
                    cleanMask(s(i):e(i)) = true;
                
               end
            
            end
                end

                
            end

            methods

            function obj = TyreModel(run1, run2)

                obj.RunData.ET = [run1.ET; run2.ET] ; % [s]
                obj.RunData.SA = [run1.SA; run2.SA]*-1 ; % [deg]
                obj.RunData.SR = [run1.SL;run2.SL]; % [-]
                obj.RunData.P = [run1.P;run2.P]*1000 ; % [Pa]
                obj.RunData.IA = [run1.IA; run2.IA] ; % [deg]
                obj.RunData.FZ = [run1.FZ;run2.FZ]; % [N]
                obj.RunData.V = [run1.V;run2.V]; % [
                obj.RunData.VX = [run1.V;run2.V];
                obj.RunData.FX = [run1.FX;run2.FX];
                obj.RunData.FY = [run1.FY;run2.FY]*-1;
                obj.RunData.MZ = [run1.MZ;run2.MZ];

            end


            function obj = createDataBins(obj, combined)

                % Pressure 
                P_tol = 0.5 ; % [psi]
                obj.P.P_8 = 0.000145038*obj.RunData.P > (8 -P_tol) & 0.000145038*obj.RunData.P < (8 +P_tol); % [0.000145038 converts Pa to psi]
                obj.P.P_10 = 0.000145038*obj.RunData.P > (10-P_tol) & 0.000145038*obj.RunData.P < (10+P_tol);
                obj.P.P_12 = 0.000145038*obj.RunData.P > (12-P_tol) & 0.000145038*obj.RunData.P < (12+P_tol);
                obj.P.P_14 = 0.000145038*obj.RunData.P > (14-P_tol) & 0.000145038*obj.RunData.P < (14+P_tol);
    
                % Inclination Angle 
                if isa(obj, 'FYclass') & combined == false 
                    IA_tol = 0.2;
                else 
                    IA_tol = 0.5;
                end

                obj.IA.IA_0 = obj.RunData.IA > - IA_tol & obj.RunData.IA < IA_tol;
                obj.IA.IA_2 = obj.RunData.IA > (2-IA_tol) & obj.RunData.IA < (2+IA_tol);
                obj.IA.IA_4 = obj.RunData.IA > (4-IA_tol) & obj.RunData.IA < (4+IA_tol);
    
                % Load
                FZ_tol = 25;
                obj.FZ.FZ_50  = 0.224809*abs(obj.RunData.FZ) > (50 -FZ_tol) & 0.224809*abs(obj.RunData.FZ) < (50 +FZ_tol);
                obj.FZ.FZ_100 = 0.224809*abs(obj.RunData.FZ) > (100-FZ_tol) & 0.224809*abs(obj.RunData.FZ) < (100+FZ_tol);
                obj.FZ.FZ_150 = 0.224809*abs(obj.RunData.FZ) > (150-FZ_tol) & 0.224809*abs(obj.RunData.FZ) < (150+FZ_tol);
                obj.FZ.FZ_200 = 0.224809*abs(obj.RunData.FZ) > (200-FZ_tol) & 0.224809*abs(obj.RunData.FZ) < (200+FZ_tol);
                obj.FZ.FZ_250 = 0.224809*abs(obj.RunData.FZ) > (250-FZ_tol) & 0.224809*abs(obj.RunData.FZ) < (250+FZ_tol);
    
                % Velocity 
                V_tol = 2 ;
                obj.V.V_15 = 0.621371*obj.RunData.V > (15-V_tol) & 0.621371*obj.RunData.V  < (15+V_tol);
                obj.V.V_25 = 0.621371*obj.RunData.V  > (25-V_tol) & 0.621371*obj.RunData.V  < (25+V_tol);
                obj.V.V_45 = 0.621371*obj.RunData.V  > (45-V_tol) & 0.621371*obj.RunData.V  < (45+V_tol);
            
                % Slip Angle 
                SA_tol = 0.5 ;
                obj.SA.SA_0 = abs(obj.RunData.SA) < SA_tol;
                obj.SA.SA_3 = abs(obj.RunData.SA) > (3-SA_tol) & abs(obj.RunData.SA) < (3+SA_tol);
                obj.SA.SA_6 = abs(obj.RunData.SA) > (6-SA_tol) & abs(obj.RunData.SA) < (6+SA_tol);

            end
        
            function obj = create_sweeps(obj, minLeng, combined)
            pressures = {'P_8', 'P_10', 'P_12', 'P_14'};
            inclinations = {'IA_0', 'IA_2', 'IA_4'};
            loads = {'FZ_50', 'FZ_100', 'FZ_150', 'FZ_200', 'FZ_250'};
        
            % Determine existing slip angles
            if combined == true
                slips = {'SA_3','SA_6'};
            else
                slips = {''}; % No slip data, use empty string
            end
        
            if nargin < 2
                minLeng = 300;
            end
        
            obj.Sweeps.sweeps = {};
            obj.Sweeps.sweep_names = {};
            obj.Sweeps.sweep_count = 0;
        
            for p_idx = 1:length(pressures)
                for ia_idx = 1:length(inclinations)
                    for sa_idx = 1:length(slips)
                        for fz_idx = 1:length(loads)
                            % Handle slip angle mask
                            if combined == true
                                sa_mask = obj.SA.(slips{sa_idx});
                                sa_val = str2double(slips{sa_idx}(4:end));
                            else
                                if isa(obj, 'FXclass')
                                    sa_mask = obj.SA.SA_0; % if no SA data, use SA = 0
                                else 
                                    sa_mask = true(size(obj.SA.SA_0)) ; 
                                end
                                
                                sa_val = NaN; 
                            end
        
                            % Sweep condition
                            sweep_condition = obj.P.(pressures{p_idx}) & ...
                                              obj.IA.(inclinations{ia_idx}) & ...
                                              sa_mask & ...
                                              obj.FZ.(loads{fz_idx}) & ...
                                              obj.V.V_25;
        
                            sweep_data = TyreModel.filterMask(sweep_condition, minLeng);
        
                            if nnz(sweep_data) > minLeng
                                obj.Sweeps.sweep_count = obj.Sweeps.sweep_count + 1;
                                obj.Sweeps.sweeps{obj.Sweeps.sweep_count} = sweep_data;
        
                                % Descriptive name
                                pressure_val = str2double(pressures{p_idx}(3:end));
                                ia_val = str2double(inclinations{ia_idx}(4:end));
                                fz_val = str2double(loads{fz_idx}(4:end));
        
                                if isnan(sa_val)
                                    obj.Sweeps.sweep_names{obj.Sweeps.sweep_count} = sprintf( ...
                                        'P=%d, IA=%d, FZ=%d (no SA data)', ...
                                        pressure_val, ia_val, fz_val);
                                else
                                    obj.Sweeps.sweep_names{obj.Sweeps.sweep_count} = sprintf( ...
                                        'P=%d, IA=%d, SA=%d, FZ=%d', ...
                                        pressure_val, ia_val, sa_val, fz_val);
                                end
        
                                fprintf('Sweep %d: %s, %d data points\n', ...
                                    obj.Sweeps.sweep_count, ...
                                    obj.Sweeps.sweep_names{obj.Sweeps.sweep_count}, ...
                                    nnz(sweep_data));
                            end
                        end
                    end
                end
            end
        end

            function obj = fitCubicSpline(obj)
                obj.splinedata = struct();
                obj.spline_fits = cell(obj.Sweeps.sweep_count);

                for i = 1:obj.Sweeps.sweep_count
                    if nnz(obj.Sweeps.sweeps{i}) > 100
                    
                        FZ_sweep = abs(obj.RunData.FZ(obj.Sweeps.sweeps{i}));
                        P_sweep = obj.RunData.P(obj.Sweeps.sweeps{i});
                        IA_sweep = obj.RunData.IA(obj.Sweeps.sweeps{i});
                        VX_sweep = obj.RunData.VX(obj.Sweeps.sweeps{i}); % Added for Mz0 function 
    
                        if isa(obj, 'FYclass')
                            slip_sweep = obj.RunData.SA(obj.Sweeps.sweeps{i});
                            force_sweep = obj.RunData.FY(obj.Sweeps.sweeps{i});
                        elseif isa(obj, 'FXclass')
                            force_sweep = obj.RunData.FX(obj.Sweeps.sweeps{i});
                            slip_sweep = obj.RunData.SR(obj.Sweeps.sweeps{i});
                        elseif isa(obj, 'MZclass')
                            force_sweep = obj.RunData.MZ(obj.Sweeps.sweeps{i});
                            slip_sweep = obj.RunData.SA(obj.Sweeps.sweeps{i});
                        end
    
                        % Sort by slips for spline fitting 
                        [slip_sorted, sort_idx] = sort(slip_sweep);
                        force_sorted = force_sweep(sort_idx);
    
                        % fit cubic smoothing spline 
                        if isa(obj, 'FXclass')
                            smoothing_param = 0.9999;
                        else
                            smoothing_param = 0.1;
                        end
    
                        obj.spline_fits{i} =csaps(slip_sorted, force_sorted, smoothing_param);
    
                        % store sweep data
                        obj.splinedata(i).slip_range = [min(slip_sorted), max(slip_sorted)];
                        obj.splinedata(i).FZ = mean(FZ_sweep);
                        obj.splinedata(i).P = mean(P_sweep);
                        obj.splinedata(i).IA = mean(IA_sweep);
                        obj.splinedata(i).name = obj.Sweeps.sweep_names{i};
                        obj.splinedata(i).VX = mean(VX_sweep);
    
                        if isa(obj, 'FYclass') | isa(obj, 'MZclass')
                            fprintf('Sweep %d (%s): SA range [%.1f, %.1f] deg, FZ = %.1fN, p = %.1fkPa, IA=%.1fdeg \n', i, obj.Sweeps.sweep_names{i}, min(slip_sorted), max(slip_sorted), mean(FZ_sweep), mean(P_sweep), mean(IA_sweep));
                        elseif isa(obj, 'FXclass')
                            fprintf('Sweep %d (%s): SR range [%.1f, %.1f] deg, FZ = %.1fN, p =%.1fkPa, IA = %.1f deg \n ', i, obj.Sweeps.sweep_names{i}, min(slip_sorted), max(slip_sorted), mean(FZ_sweep), mean(P_sweep), mean(IA_sweep));
                        end
                    end
                    end
            end

            function obj = fitCombinedCubicSpline(obj)
                obj.spline_fits = cell(1, obj.Sweeps.sweep_count);
                obj.splinedata = struct();
                
                for i = 1:obj.Sweeps.sweep_count 
                    if nnz(obj.Sweeps.sweeps{i}) > 100

                        % extract data for this sweep 
                        SR_sweep = obj.RunData.SR(obj.Sweeps.sweeps{i});
                        SA_sweep = obj.RunData.SA(obj.Sweeps.sweeps{i});
                        FZ_sweep = abs(obj.RunData.FZ(obj.Sweeps.sweeps{i}));
                        P_sweep = obj.RunData.P(obj.Sweeps.sweeps{i});
                        IA_sweep = obj.RunData.IA(obj.Sweeps.sweeps{i});
                        VX_sweep = obj.RunData.VX(obj.Sweeps.sweeps{i});

                        if isa(obj, 'FYclass')
                            force_sweep = obj.RunData.FY(obj.Sweeps.sweeps{i});
                        elseif isa(obj, 'FXclass')
                            force_sweep = obj.RunData.FX(obj.Sweeps.sweeps{i});
                        elseif isa(obj, 'MZclass')
                            force_sweep = obj.RunData.MZ(obj.Sweeps.sweeps{i}) ; 
                        end

                        [SR_sorted, sort_idx] = sort(SR_sweep);
                        force_sorted = force_sweep(sort_idx);

                        [SR_unique, unique_idx] = unique(SR_sorted);
                        F_unique = force_sorted(unique_idx);

                        smoothing_param = 0.9999;
                        obj.spline_fits{i} = csaps(SR_unique, F_unique, smoothing_param) ;

                        % store sweep data
                        obj.splinedata(i).SR_range = [min(SR_unique), max(SR_unique)];
                        obj.splinedata(i).SA = mean(SA_sweep);
                        obj.splinedata(i).FZ = mean(FZ_sweep);
                        obj.splinedata(i).P = mean(P_sweep);
                        obj.splinedata(i).IA = mean(IA_sweep);
                        obj.splinedata(i).name = obj.Sweeps.sweep_names{i};
                        obj.splinedata(i).VX = mean(VX_sweep);

                        if isa(obj, 'MZclass')
                            fprintf('Sweep %d (%s): SR range [%.1f, %.1f]deg, FZ=%.1fN, p=%.1fkPa, IA=%.1fdeg \n, SA=%.1fdeg\n, VX=%.1fN\n', i, obj.Sweeps.sweep_names{i}, min(SR_unique), max(SR_unique), mean(FZ_sweep), mean(P_sweep), mean(IA_sweep), mean(SA_sweep), mean(VX_sweep))
                        else
                            fprintf('Sweep %d (%s): SR range [%.1f, %.1f]deg, FZ=%.1fN, p=%.1fkPa, IA=%.1fdeg \n, SA=%.1fdeg\n', i, obj.Sweeps.sweep_names{i}, min(SR_unique), max(SR_unique), mean(FZ_sweep), mean(P_sweep), mean(IA_sweep), mean(SA_sweep))
                    
                        end
                    end
                end
            end

            function obj = generateSplinePoints(obj) 
                n = 100 ; % number of points per spline 
                obj.evaluation_points.all_slip = [] ;
                obj.evaluation_points.all_Fz = [] ;
                obj.evaluation_points.all_pressure = [] ;
                obj.evaluation_points.all_inclangle = [] ;
                obj.evaluation_points.all_force_spline = [];
                obj.evaluation_points.sweep_indices = [];
                obj.evaluation_points.all_Vcx = []; % Added for Mz0 function 

                for i = 1:obj.Sweeps.sweep_count 
                    if ~isempty(obj.splinedata(i).slip_range)
                        slip_eval = linspace(obj.splinedata(i).slip_range(1), obj.splinedata(i).slip_range(2), n) ; 
                        
                        force_spline_eval = fnval(obj.spline_fits{i}, slip_eval);

                        % store the evaluation points 
                        
                        if isa(obj, "FXclass")
                            obj.evaluation_points.all_slip = [obj.evaluation_points.all_slip; slip_eval'] ;
                        else
                            obj.evaluation_points.all_slip = [obj.evaluation_points.all_slip; deg2rad(slip_eval')] ;
                        end
                        obj.evaluation_points.all_Fz = [obj.evaluation_points.all_Fz; obj.splinedata(i).FZ * ones(n,1)] ;
                        obj.evaluation_points.all_pressure = [obj.evaluation_points.all_pressure ; obj.splinedata(i).P * ones(n,1)] ;
                        obj.evaluation_points.all_inclangle = [obj.evaluation_points.all_inclangle ; deg2rad(obj.splinedata(i).IA)* ones(n,1)] ;
                        obj.evaluation_points.all_force_spline = [obj.evaluation_points.all_force_spline; force_spline_eval'];
                        obj.evaluation_points.sweep_indices = [obj.evaluation_points.sweep_indices; i*ones(n, 1) ];
                        obj.evaluation_points.all_Vcx = [obj.evaluation_points.all_Vcx; obj.splinedata(i).VX * ones(n,1)];

                    end
                end

                fprintf('Total evaluation points for optimisation: %d\n', length(obj.evaluation_points.all_force_spline))

            end

            function obj = generateCombinedSplinePoints(obj)
                n = 100; % number of points per spline 
                obj.evaluation_points.all_longslip = [];
                obj.evaluation_points.all_slipangl = [];
                obj.evaluation_points.all_Fz = [];
                obj.evaluation_points.all_pressure = [];
                obj.evaluation_points.all_inclangl = [];
                obj.evaluation_points.all_force_spline = [];
                obj.evaluation_points.sweep_indices = [];
                obj.evaluation_points.all_Vcx = [];

                for i = 1:obj.Sweeps.sweep_count
                    if ~isempty(obj.splinedata(i).SR_range)
                        % create evenly spaced slip angles across the sweep
                        % range 
                        SR_eval = linspace(obj.splinedata(i).SR_range(1), obj.splinedata(i).SR_range(2), n);

                        force_spline_eval = fnval(obj.spline_fits{i}, SR_eval);

                        % store the evaluation points
                        obj.evaluation_points.all_longslip = [obj.evaluation_points.all_longslip; SR_eval'];
                        obj.evaluation_points.all_slipangl = [obj.evaluation_points.all_slipangl; deg2rad(obj.splinedata(i).SA) * ones(n,1) ];
                        obj.evaluation_points.all_Fz = [obj.evaluation_points.all_Fz; obj.splinedata(i).FZ * ones(n,1)] ;
                        obj.evaluation_points.all_pressure = [obj.evaluation_points.all_pressure ; obj.splinedata(i).P * ones(n,1)] ;
                        obj.evaluation_points.all_inclangl = [obj.evaluation_points.all_inclangl ; deg2rad(obj.splinedata(i).IA)* ones(n,1)] ;
                        obj.evaluation_points.all_force_spline = [obj.evaluation_points.all_force_spline; force_spline_eval'];
                        obj.evaluation_points.sweep_indices = [obj.evaluation_points.sweep_indices; i*ones(n, 1) ];
                        obj.evaluation_points.all_Vcx = [obj.evaluation_points.all_Vcx; obj.splinedata(i).VX * ones(n,1)];
                    end
                end

                fprintf('Total evaluation points for optimisation: %d\n', length(obj.evaluation_points.all_force_spline))
            end

            function [obj, updated_p] = magicFormulaParamterOptimisation(obj,p)
                updated_p = p;

                 obj = obj.optimisationParameters() ;
                 obj = obj.initialGuess();
                 
                 % find the initial fit of the magic formula 
                 if isa(obj, 'FYclass')
                    
                    [F0_initial, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                        FYclass.Fy0(obj.p, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle);
                    
                    residual_fun = @(x) FYclass.residual_fy0(x, obj.p, obj.params.params_to_optimize, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle, obj.evaluation_points.all_force_spline, obj.evaluation_points.sweep_indices);

                 elseif isa(obj, 'FXclass')
                    
                    [F0_initial, ~, ~, ~, ~, ~, ~] = FXclass.Fx0(obj.p, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle);
                    
                    residual_fun = @(x) FXclass.residual_fx0(x, obj.p, obj.params.params_to_optimize, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle, obj.evaluation_points.all_force_spline, obj.evaluation_points.sweep_indices);

                 elseif isa(obj, 'MZclass')
                    [F0_initial, ~, ~, ~, ~, ~, ~, ~, ~] = MZclass.Mz0(obj.p, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle, obj.evaluation_points.all_Vcx);
                    
                    residual_fun = @(x) MZclass.residual_mz0(x, obj.p, obj.params.params_to_optimize, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle, obj.evaluation_points.all_Vcx, obj.evaluation_points.all_force_spline, obj.evaluation_points.sweep_indices);
                 end

                 % Run Optimisation 
                 fprintf('Starting optimisation with %d spline curves...\n', obj.Sweeps.sweep_count);
                 [x_opt, obj.evaluation_points.resnorm, residual, obj.evaluation_points.exitflag, output] = ...
                     lsqnonlin(residual_fun, obj.params.x0, obj.params.lb, obj.params.ub, obj.params.options);

                     for i = 1:length(obj.params.params_to_optimize)

                        % Update parameters to optimise 
                        param_name = obj.params.params_to_optimize{i};
                        obj.p.(param_name) = x_opt(i);
                        updated_p.(param_name) = x_opt(i);
                     end

                if isa(obj, 'FYclass')
                % Calculate Magic Formula fit for Fy0
                    [obj.magic_formula.F0_fit, obj.magic_formula.mu0_fit, obj.magic_formula.dfz, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                    FYclass.Fy0(obj.p, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle);  
                elseif isa(obj, 'FXclass')
                % Calculate Magic Formula fit with optimized parameters
                    [obj.magic_formula.F0_fit, obj.magic_formula.mu0_fit, Cx_fit, Dx_fit, Ex_fit, dfz_fit, Kxk_fit] = ...
                    FXclass.Fx0(obj.p, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle);
        
                elseif isa(obj, 'MZclass')
                % Calculate Magic Formula fit with optimized parameters
                    [obj.magic_formula.F0_fit, Bt_fit, Ct_fit, Dt_fit, Et_fit, Br_fit, Cr_fit, Dr_fit, alphaCos_fit, R0_fit, alphar_fit, alphat_fit] = ...
                    MZclass.Mz0(obj.p, obj.evaluation_points.all_slip, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangle, obj.evaluation_points.all_Vcx);
                   
                end
              

            end
       

            function [obj, updated_p] = magicFormulaCombinedParameterOptimisation(obj, p) 
                updated_p = p;
                obj = obj.combinedoptimisationParameters(); % creates obj.p for all coefficients related to Fy 

                obj = obj.initialGuess() ;
                parameterfields = fieldnames(p);

                for i = 1:numel(parameterfields) % updates all the coefficients that have been optimised previously (single slip) 
                    field = parameterfields{i};
                    if isfield(obj.p, field)
                        obj.p.(field) = p.(field);
                    end
                end

                % find the initial fit to the spline curve 
                if isa(obj, 'FYclass')
                    [F_initial, ~, ~, ~, ~, ~, ~] = FYclass.Fy(obj.p, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl);
                    residual_fun = @(x) FYclass.residual_fy(x, obj.p, obj.params.params_to_optimize, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl, obj.evaluation_points.all_force_spline, obj.evaluation_points.sweep_indices);
                elseif isa(obj, 'FXclass')
                    [F_initial, ~, ~, ~, ~] = FXclass.Fx(obj.p, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl);
                    residual_fun = @(x) FXclass.residual_fx(x, obj.p, obj.params.params_to_optimize, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl, obj.evaluation_points.all_force_spline, obj.evaluation_points.sweep_indices);
                
                elseif isa(obj, 'MZclass')
                    F_initial = MZclass.Mz(obj.p, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl, obj.evaluation_points.all_Vcx, obj.evaluation_points.all_Fx, obj.evaluation_points.all_Fy);
                    residual_fun = @(x) MZclass.residual_mz(x, obj.p, obj.params.params_to_optimize, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl, obj.evaluation_points.all_Vcx, obj.evaluation_points.all_Fx, obj.evaluation_points.all_Fy, obj.evaluation_points.all_force_spline, obj.evaluation_points.sweep_indices);
                end


                % Run optimisation 
                fprintf('Starting combined optimisation with %d spline curves... \n', obj.Sweeps.sweep_count);
                [x_opt, obj.evaluation_points.resnorm, residual, obj.evaluation_points.exitflag, output] = ...
                    lsqnonlin(residual_fun, obj.params.x0, obj.params.lb, obj.params.ub, obj.params.options);
                
                for i = 1:length(obj.params.params_to_optimize)
                    param_name = obj.params.params_to_optimize{i};
                    obj.p.(param_name) = x_opt(i);
                    updated_p.(param_name) = x_opt(i);
                end

                % Calculate the Magic Formula Fit 
                if isa(obj, 'FYclass')
                    [obj.magic_formula.F_fit, ~, ~, ~, ~, ~, ~] = FYclass.Fy(obj.p, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl);
                elseif isa(obj, 'FXclass')
                    [obj.magic_formula.F_fit, ~, ~, ~, ~] = FXclass.Fx(obj.p, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl);
                elseif isa(obj, 'MZclass')
                    obj.magic_formula.F_fit = MZclass.Mz(obj.p, obj.evaluation_points.all_longslip, obj.evaluation_points.all_slipangl, obj.evaluation_points.all_Fz, obj.evaluation_points.all_pressure, obj.evaluation_points.all_inclangl, obj.evaluation_points.all_Vcx, obj.evaluation_points.all_Fx, obj.evaluation_points.all_Fy);
                end
               
            end

        % function obj = plotting_pressure_camber(obj)
        %     load_levels = [50, 100, 150, 200, 250];
        %     for load_idx = 1:length(load_levels)
        %         current_load = 4.448*load_levels(load_idx);
        % 
        %         Find all sweeps with this load
        %         load_sweeps = [];
        %         for i = 1:obj.Sweeps.sweep_count
        %             if abs(obj.splinedata(i).FZ - current_load) < 30
        %                 load_sweeps = [load_sweeps, i];
        %             end
        %         end
        % 
        %         if ~isempty(load_sweeps)
        %             Create figure for this load level
        %             figure('Name', sprintf('FZ=%dN Sweeps', current_load));
        %             legend_handles = gobjects(3,1);
        %             legend_created = false;
        % 
        %             Calculate subplot layout
        %             num_sweeps = length(load_sweeps);
        %             rows = ceil(sqrt(num_sweeps));
        %             cols = ceil(num_sweeps / rows);
        % 
        %             for plot_idx = 1:num_sweeps
        %                 i = load_sweeps(plot_idx);
        %                 subplot(rows, cols, plot_idx);
        %                 hold on;
        % 
        %                 Get evaluation points for this sweep
        %                 sweep_mask = (obj.evaluation_points.sweep_indices == i);
        %                 if any(sweep_mask)
        %                     if isa(obj, "FYclass")
        %                         slip_plot = rad2deg(obj.evaluation_points.all_slip(sweep_mask));
        %                     elseif isa(obj, "FXclass")
        %                         slip_plot = obj.evaluation_points.all_slip(sweep_mask);
        %                     elseif isa(obj, "MZclass")
        %                         slip_plot = rad2deg(obj.evaluation_points.all_slip(sweep_mask));
        %                     end
        % 
        %                     Plot spline curve
        %                     h1 = plot(slip_plot, obj.evaluation_points.all_force_spline(sweep_mask), 'g-', 'LineWidth', 2, 'DisplayName', 'Spline');
        % 
        %                     Plot Magic Formula fit
        %                     h2 = plot(slip_plot, obj.magic_formula.F0_fit(sweep_mask), 'r--', 'LineWidth', 2, 'DisplayName', 'Magic Formula');
        % 
        %                     Plot raw data (if available)
        %                     h3 = [];
        %                     if i <= length(obj.Sweeps.sweeps)
        %                         if isa(obj, "FYclass")
        %                             raw_slip = obj.RunData.SA(obj.Sweeps.sweeps{i});
        %                             raw_force = obj.RunData.FY(obj.Sweeps.sweeps{i});
        %                         elseif isa(obj, "FXclass")
        %                             raw_slip = obj.RunData.SR(obj.Sweeps.sweeps{i});
        %                             raw_force = obj.RunData.FX(obj.Sweeps.sweeps{i});
        % 
        %                         elseif isa(obj, "MZclass")
        %                             raw_slip = obj.RunData.SA(obj.Sweeps.sweeps{i});
        %                             raw_force = obj.RunData.MZ(obj.Sweeps.sweeps{i});
        %                         end
        %                         h3 = plot(raw_slip, raw_force, 'b.', 'MarkerSize', 2, 'DisplayName', 'Raw Data');
        %                     end
        % 
        %                     if ~legend_created
        %                         legend_handles(1) = h1;
        %                         legend_handles(2) = h2;
        %                         if ~isempty(h3)
        %                             legend_handles(3) = h3;
        %                         end
        %                         legend_created = true;
        %                     end
        % 
        % 
        %                     Create title with all parameters
        %                     title_str = sprintf('P=%d, IA=%d, FZ=%d', ...
        %                         round(obj.splinedata(i).P), round(obj.splinedata(i).IA), round(obj.splinedata(i).FZ));
        %                     title(title_str, 'FontSize', 10);
        % 
        % 
        %                     if isa(obj, "FYclass")
        %                         ylabel('F_y (N)');
        %                         xlabel('Slip Angle (deg)');
        %                     elseif isa(obj, "FXclass")
        %                         ylabel('F_x (N)');
        %                         xlabel('Slip Ratio (-)');
        %                     elseif isa(obj, "MZclass")
        %                         ylabel('M_z (N)');
        %                         xlabel('Slip Angle (deg)');
        %                     end
        % 
        %                     grid on;
        % 
        %                     hold off
        %                 end
        %             end
        % 
        %             Add overall title
        %             sgtitle(sprintf('FZ = %dN - Pressure and Camber Variations', current_load), 'FontSize', 14, 'FontWeight', 'bold');
        % 
        %             figure_legend = legend(legend_handles, {'Spline','Magic Formula','Raw Data'}, ...
        %                            'Orientation','horizontal', 'Location','southoutside');
        % 
        %             pos = figure_legend.Position;   % current normalized position
        %             pos(2) = 0.01;                  % move legend close to bottom
        %             pos(1) = 0.25;                  % center horizontally
        %             pos(3) = 0.5;                   % width
        %             figure_legend.Position = pos;
        %         end
        %     end
        % end 
        function draw_pressure_camber(obj,current_load)
        
            % Find all sweeps with this load
            load_sweeps = [];
            for i = 1:obj.Sweeps.sweep_count
                if abs(obj.splinedata(i).FZ - current_load) < 30
                    load_sweeps = [load_sweeps, i];
                end
            end
        
            if isempty(load_sweeps)
                warning('No sweeps found for this load')
                return
            end
        
            figure('Name', sprintf('FZ=%dN Sweeps', current_load));
            legend_handles = gobjects(3,1);
            legend_created = false;
        
            % subplot layout
            num_sweeps = length(load_sweeps);
            rows = ceil(sqrt(num_sweeps));
            cols = ceil(num_sweeps / rows);
        
            for plot_idx = 1:num_sweeps
                i = load_sweeps(plot_idx);
                subplot(rows, cols, plot_idx);
                hold on;
        
                sweep_mask = (obj.evaluation_points.sweep_indices == i);
        
                if any(sweep_mask)
        
                    if isa(obj,"FYclass") || isa(obj,"MZclass")
                        slip_plot = rad2deg(obj.evaluation_points.all_slip(sweep_mask));
                    else
                        slip_plot = obj.evaluation_points.all_slip(sweep_mask);
                    end
        
                    h1 = plot(slip_plot, obj.evaluation_points.all_force_spline(sweep_mask),'g-','LineWidth',2);
                    h2 = plot(slip_plot, obj.magic_formula.F0_fit(sweep_mask),'r--','LineWidth',2);
        
                    h3 = [];
        
                    if i <= length(obj.Sweeps.sweeps)
        
                        if isa(obj,"FYclass")
                            raw_slip = obj.RunData.SA(obj.Sweeps.sweeps{i});
                            raw_force = obj.RunData.FY(obj.Sweeps.sweeps{i});
        
                        elseif isa(obj,"FXclass")
                            raw_slip = obj.RunData.SR(obj.Sweeps.sweeps{i});
                            raw_force = obj.RunData.FX(obj.Sweeps.sweeps{i});
        
                        elseif isa(obj,"MZclass")
                            raw_slip = obj.RunData.SA(obj.Sweeps.sweeps{i});
                            raw_force = obj.RunData.MZ(obj.Sweeps.sweeps{i});
                        end
        
                        h3 = plot(raw_slip, raw_force,'b.','MarkerSize',2);
                    end
        
                    if ~legend_created
                        legend_handles(1)=h1;
                        legend_handles(2)=h2;
                        if ~isempty(h3); legend_handles(3)=h3; end
                        legend_created=true;
                    end
        
                    title(sprintf('P=%d, IA=%d, FZ=%d',...
                        round(obj.splinedata(i).P),...
                        round(obj.splinedata(i).IA),...
                        round(obj.splinedata(i).FZ)))
        
                    if isa(obj,"FYclass")
                        ylabel('F_y (N)'); xlabel('Slip Angle (deg)');
                    elseif isa(obj,"FXclass")
                        ylabel('F_x (N)'); xlabel('Slip Ratio (-)');
                    else
                        ylabel('M_z (Nm)'); xlabel('Slip Angle (deg)');
                    end
        
                    grid on
                    hold off
                end
            end
        
            sgtitle(sprintf('FZ = %dN - Pressure and Camber Variations',current_load),'FontSize',14,'FontWeight','bold');
        
            legend(legend_handles,{'Spline','Magic Formula','Raw Data'},...
                'Orientation','horizontal','Location','southoutside');
        end

    function obj = buildPlotLibrary(obj,varargin)
        % Usage:
        %   obj = obj.buildPlotLibrary("combined","yes")
        %   obj = obj.buildPlotLibrary("combined","no")
        
            useCombined = false;
        
            if nargin > 1
                for k = 1:2:length(varargin)
                    if strcmpi(varargin{k},"combined")
                        if strcmpi(varargin{k+1},"yes")
                            useCombined = true;
                        end
                    end
                end
            end
            obj.StoredPlots = struct();
            obj.StoredPlots.plots = {};
            obj.StoredPlots.loads = [];

            load_levels = [50,100,150,200,250];
        
            pc = 1;
            cc = 1;
        
            for load_idx = 1:length(load_levels)
        
                current_load = 4.44822 * load_levels(load_idx);
        
                % Detect if load exists
                load_exists = false;
                for i = 1:obj.Sweeps.sweep_count
                    if ~isempty(obj.splinedata(i).FZ) && abs(obj.splinedata(i).FZ - current_load) < 30
                        load_exists = true;
                        break
                    end
                end
        
                if ~load_exists
                    continue
                end
  
                if useCombined
                    obj.StoredPlots.plots{cc} = ...
                        @() obj.draw_combined_plots(current_load);
        
                    obj.StoredPlots.loads(cc) = current_load;
                    cc = cc + 1;
                else % for non combined
                    obj.StoredPlots.plots{pc} = ...
                    @() obj.draw_pressure_camber(current_load);
        
                    obj.StoredPlots.loads(pc) = current_load;
                    pc = pc + 1;

                end
            end
        
        end
        
        function obj = draw_combined_plots(obj, current_load)

            % Find sweeps with this load
            load_sweeps = [];
            for i = 1:obj.Sweeps.sweep_count
                if ~isempty(obj.splinedata(i).FZ) && abs(obj.splinedata(i).FZ - current_load) < 30
                    load_sweeps = [load_sweeps, i];
                end
            end
        
            if isempty(load_sweeps)
                warning('No sweeps found for this load')
                return
            end
        
            figure('Name', sprintf('FZ=%dN Sweeps', round(current_load)),...
                   'Position',[100 100 1200 800]);
        
            legend_handles = gobjects(3,1);
            legend_created = false;
        
            num_sweeps = length(load_sweeps);
            rows = ceil(sqrt(num_sweeps));
            cols = ceil(num_sweeps / rows);
        
            for plot_idx = 1:num_sweeps
                i = load_sweeps(plot_idx);
                subplot(rows,cols,plot_idx);
                hold on;
        
                sweep_mask = (obj.evaluation_points.sweep_indices == i);
        
                if any(sweep_mask)
        
                    SR = obj.evaluation_points.all_longslip(sweep_mask);
        
                    h1 = plot(SR,obj.evaluation_points.all_force_spline(sweep_mask),'g-','LineWidth',2);
                    h2 = plot(SR,obj.magic_formula.F_fit(sweep_mask),'r--','LineWidth',2);
        
                    h3 = [];
                    if i <= length(obj.Sweeps.sweeps) && ~isempty(obj.Sweeps.sweeps{i})
        
                        raw_SR = obj.RunData.SR(obj.Sweeps.sweeps{i});
        
                        if isa(obj,'FYclass')
                            raw_F = obj.RunData.FY(obj.Sweeps.sweeps{i});
                        elseif isa(obj,'FXclass')
                            raw_F = obj.RunData.FX(obj.Sweeps.sweeps{i});
                        else
                            raw_F = obj.RunData.MZ(obj.Sweeps.sweeps{i});
                        end
        
                        h3 = plot(raw_SR,raw_F,'b.','MarkerSize',2);
                    end
        
                    if ~legend_created
                        legend_handles(1)=h1;
                        legend_handles(2)=h2;
                        if ~isempty(h3); legend_handles(3)=h3; end
                        legend_created=true;
                    end
        
                    title(sprintf('P=%.1fkPa, IA=%.1fdeg, SA=%.1fdeg, FZ=%.1fN',...
                        obj.splinedata(i).P,...
                        obj.splinedata(i).IA,...
                        obj.splinedata(i).SA,...
                        obj.splinedata(i).FZ),'FontSize',10);
        
                    xlabel('Slip Ratio');
        
                    if isa(obj,'FXclass')
                        ylabel('F_x (N)');
                    elseif isa(obj,'FYclass')
                        ylabel('F_y (N)');
                    else
                        ylabel('M_z (Nm)');
                    end
        
                    grid on
                    hold off
                end
            end
            % Create legend at very bottom
            lgd = legend(legend_handles, {'Spline','Magic Formula','Raw Data'}, 'Orientation', 'horizontal');
            lgd.Units = 'normalized';
            lgd.Position = [0.25, 0.01, 0.5, 0.03];  % [x, y, width, height]
                        
        end

        function showPlot(obj,load_value)
            if ~isfield(obj.StoredPlots,'plots')
                error('Run BuildPlotLibrary first')
            end

            [~,idx] = min(abs(obj.StoredPlots.loads - load_value));
        
            obj.StoredPlots.plots{idx}();
        
        end
        
           
       %  function obj = plotting_for_combined(obj)
       %      load_levels = [50, 100, 150, 200, 250];
       % 
       %      for load_idx = 1:length(load_levels)
       %          current_load = 4.44822 * load_levels(load_idx); % convert Lbs to N
       % 
       %          % Find all sweeps with this load 
       %          load_sweeps = []; 
       %          for i = 1:obj.Sweeps.sweep_count
       %              if ~isempty(obj.splinedata(i).FZ) && abs(obj.splinedata(i).FZ - current_load) < 30 
       %                  load_sweeps = [load_sweeps, i];
       %              end
       %          end
       % 
       %          if ~isempty(load_sweeps) 
       %              % Create figure for this load level 
       %              figure('Name', sprintf('FZ=%dN Sweeps', round(current_load)), 'Position', [100, 100, 1200, 800]);
       % 
       %              legend_handles = gobjects(3,1);
       %              legend_created = false;
       % 
       %              % Calculate subplot layout
       %              num_sweeps = length(load_sweeps);
       %              rows = ceil(sqrt(num_sweeps));
       %              cols = ceil(num_sweeps / rows);
       % 
       %              for plot_idx = 1:num_sweeps
       %                  i = load_sweeps(plot_idx);
       %                  subplot(rows, cols, plot_idx);
       %                  hold on;
       % 
       %                  % Get evaluation points for this sweep
       %                  sweep_mask = (obj.evaluation_points.sweep_indices == i);
       %                  if any(sweep_mask)
       %                      SR = obj.evaluation_points.all_longslip(sweep_mask);
       % 
       %                      % Plot spline curve
       %                      h1 = plot(SR, obj.evaluation_points.all_force_spline(sweep_mask), 'g-', 'LineWidth', 2, 'DisplayName', 'Spline');
       %                      % Plot Magic Formula fit
       %                      h2 = plot(SR, obj.magic_formula.F_fit(sweep_mask), 'r--', 'LineWidth', 2, 'DisplayName', 'Magic Formula');
       % 
       %                      % Plot raw data (if available)
       %                      h3 = [];
       %                      if i <= length(obj.Sweeps.sweeps) && ~isempty(obj.Sweeps.sweeps{i})
       %                          raw_SR = obj.RunData.SR(obj.Sweeps.sweeps{i});
       %                          if isa(obj, 'FYclass')
       %                              raw_F = obj.RunData.FY(obj.Sweeps.sweeps{i});
       %                          elseif isa(obj, 'FXclass')
       %                              raw_F = obj.RunData.FX(obj.Sweeps.sweeps{i});
       %                          elseif isa(obj, 'MZclass')
       %                              raw_F = obj.RunData.MZ(obj.Sweeps.sweeps{i});
       %                          end
       % 
       %                          h3 = plot(raw_SR, raw_F, 'b.', 'MarkerSize', 2, 'DisplayName', 'Raw Data');
       %                      end
       % 
       %                      if ~legend_created
       %                          legend_handles(1) = h1;
       %                          legend_handles(2) = h2;
       %                          if ~isempty(h3)
       %                              legend_handles(3) = h3;
       %                          end
       %                          legend_created = true;
       %                      end
       % 
       %                      % Create title with all parameters 
       %                      title_str = sprintf('P=%.1fkPa, IA=%.1fdeg, SA=%.1fdeg, FZ=%.1fN', ...
       %                          obj.splinedata(i).P, obj.splinedata(i).IA, obj.splinedata(i).SA, obj.splinedata(i).FZ);
       %                      title(title_str, 'FontSize', 10);
       % 
       % 
       %                      xlabel('Slip Ratio');
       %                      if isa(obj, 'FXclass')
       %                          ylabel('F_x (N)');
       %                      elseif isa(obj, 'FYclass')
       %                          ylabel('F_y (N)');
       %                      end
       %                      grid on;
       % 
       % 
       %                      hold off;
       %                  end
       %              end
       %              figure_legend = legend(legend_handles, {'Spline','Magic Formula','Raw Data'}, ...
       % 'Orientation','horizontal', 'Location','southoutside');
       %                                  pos = figure_legend.Position;   % current normalized position
       %              pos(2) = 0.01;                  % move legend close to bottom
       %              pos(1) = 0.25;                  % center horizontally
       %              pos(3) = 0.5;                   % width
       %              figure_legend.Position = pos;
       %          end
       %      end
       %  end

        function obj = resultsAnalysis(obj, combined)
            fprintf('\n=== OPTIMISATION RESULTS ===\n');
            fprintf('Total Sweeps: %d\n', obj.Sweeps.sweep_count);
            fprintf('Total Evalulation Points: %d\n', length(obj.evaluation_points.all_force_spline));
            fprintf('Exit Flag: %d\n', obj.evaluation_points.exitflag);
            fprintf('Resiudal Norm: %d\n', obj.evaluation_points.resnorm);

            % Finding R squared of the residuals 
            if combined == false
                fit = obj.magic_formula.F0_fit;
            else
                fit = obj.magic_formula.F_fit;

            end
            
            ss_res = sum((obj.evaluation_points.all_force_spline - fit).^2);
            ss_tot = sum((obj.evaluation_points.all_force_spline - mean(obj.evaluation_points.all_force_spline)).^2);
            r_squared = 1 - (ss_res / ss_tot);
            fprintf('Overall R^2 of the fit: %.4f\n', r_squared);

            % Sweep by sweep performance 
            fprintf('\n=== SWEEP BY SWEEP PERFORMANCE ===\n');
            for i = 1:obj.Sweeps.sweep_count
                sweep_mask = (obj.evaluation_points.sweep_indices == i);
                if any(sweep_mask)
                    sweep_spline = obj.evaluation_points.all_force_spline(sweep_mask);
                    sweep_fit = fit(sweep_mask);
                    ss_res = sum((sweep_spline - sweep_fit).^2);
                    ss_tot = sum((sweep_spline - mean(sweep_spline)).^2);
                    r_squared = 1- (ss_res/ss_tot);
                    rmse = sqrt(mean((sweep_spline - sweep_fit).^2));

                    fprintf('Sweep %d (%s): R^2 = %.4f, RMSE=%.2fN\n', i, obj.splinedata(i).name, r_squared, rmse);
                end
            end

            fprintf('\n===OPTIMISED PARAMETERS===\n')
            for i = 1:length(obj.params.params_to_optimize)
                param_name = obj.params.params_to_optimize{i};
                fprintf('%s = %.6f\n', param_name, obj.p.(param_name));
            end
        end

        obj = plot_Fy_Fz_SA(obj)

    end
end

            



                           