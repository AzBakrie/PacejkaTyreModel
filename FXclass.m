classdef FXclass < TyreModel & Parameters

    properties
        params 
    end

    methods(Static)
        function [Fx0,mux,Cx,Dx,Ex,dfz,Kxk] = Fx0(p,longslip,Fz,pressure,inclangl)
            %MAGICFORMULA.V61.FX0
            % TODO: implement for turnslip (ZETAs):
            ZETA1 = 1;
            
            % (4.E1, 4.E2a)
            Fz0_ = p.FNOMIN.*p.LFZO;
            dfz = (Fz-Fz0_)./Fz0_;
            
            % (4.E2b) 
            pi0 = p.NOMPRES;
            dpi = (pressure-pi0)./pi0;
            
            % (4.E11) 
            Cx = p.PCX1.*p.LCX;
            
            % (4.E13) 
            % Note: Asterisk variant of LMUX is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            mux = p.LMUX.*(p.PDX1+p.PDX2.*dfz).*(1+p.PPX3.*dpi+p.PPX4.*dpi.^2).*(1-p.PDX3.*inclangl.^2);
            
            % (4.E12) 
            Dx = mux.*Fz*ZETA1;
            
            % (4.E17) 
            SHx = p.LHX.*(p.PHX1+p.PHX2.*dfz);
            
            % (4.E8)
            % Note: Asterisk variant of LMUX is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            AMU = 10;
            LMUX_ = AMU.*p.LMUX./(1+(AMU-1).*p.LMUX);
            
            % (4.E18)
            SVx = ZETA1*LMUX_.*p.LVX.*Fz.*(p.PVX1+p.PVX2.*dfz);
            
            % (4.E10) 
            kappax = longslip + SHx;
            kappaxSgn = sign(kappax);
            
            % (4.E14) 
            Ex = p.LEX.*(p.PEX1+p.PEX2.*dfz+p.PEX3.*dfz.^2).*(1-p.PEX4.*kappaxSgn);
            
            % (4.E15)
            Kxk = p.LKX.*Fz.*(p.PKX1+p.PKX2.*dfz).*exp(p.PKX3.*dfz).*(1+p.PPX1.*dpi+p.PPX2.*dpi.^2);
            
            % (4.E16) 
            Bx = Kxk./(Cx.*Dx+eps(Kxk));
            
            % (4.E9)
            Fx0 = Dx.*sin(Cx.*atan(Bx.*kappax-Ex.*(Bx.*kappax-atan(Bx.*kappax))))+SVx;

        end

        function [Fx,mux,Gxa,Bxa,Exa] = Fx(p,longslip,slipangl,Fz,pressure,inclangl)
            %MAGICFORMULA.V61.FX
            [Fx0_,mux,~,~,~,dfz] = FXclass.Fx0(p,longslip,Fz,pressure,inclangl);
            
            % (4.E3)
            alphaAst = tan(slipangl);
            
            % (4.E4)
            gammaAst = sin(inclangl);
            
            % (4.E55)
            Cxa = p.RCX1;
            
            % (4.E56)
            Exa = p.REX1 + p.REX2.*dfz;
            
            % (4.E57)
            SHxa = p.RHX1;
            
            % (4.E53)
            alphaS = alphaAst + SHxa;
            
            % (4.E54)
            Bxa = p.LXAL*(p.RBX1+p.RBX3*gammaAst.^2).*cos(atan(p.RBX2.*longslip));
            
            % (4.E52)
            Gxa0 = cos(Cxa.*atan(Bxa.*SHxa-Exa.*(Bxa.*SHxa-atan(Bxa.*SHxa))));
            
            % (4.E51)
            Gxa = cos(Cxa.*atan(Bxa.*alphaS-Exa.*(Bxa.*alphaS-atan(Bxa.*alphaS))))./Gxa0;
            
            % (4.E50)
            Fx = Gxa.*Fx0_;
        end
        
        %% Residual function with sweep normalization for Fy0
        function residuals = residual_fx0(x, p, param_names, ...
            longslip, Fz, pressure, inclangl, Fx_spline_target, sweep_indices)
        
        % Update parameter structure
        p_current = p;
        for i = 1:length(param_names)
            p_current.(param_names{i}) = x(i);
        end
        
        try
            % Calculate Magic Formula output
            [Fx0_pred, ~, ~, ~, ~, ~, ~] = ...
                FXclass.Fx0(p_current, longslip, Fz, pressure, inclangl);
            
            % Calculate base residuals
            residuals = Fx0_pred - Fx_spline_target;
           
            % Normalize by sweep to ensure equal contribution
            unique_sweeps = unique(sweep_indices);
            for sweep_idx = unique_sweeps'
                sweep_mask = (sweep_indices == sweep_idx);
                if any(sweep_mask)
                    sweep_max_force = max(abs(Fx_spline_target(sweep_mask)));
                    if sweep_max_force > 0 
                        residuals(sweep_mask) = residuals(sweep_mask)/sweep_max_force; 
                
                    end
                end

            end
           
            
        catch ME
            fprintf('Error in Fx0: %s\n', ME.message);
            residuals = 1000 * ones(size(Fx_spline_target));
        end
        end
        %% Residual function with sweep normalization
        function residuals = residual_fx(x, p, param_names, ...
            longslip, slipangl, Fz, pressure, inclangl, Fx_spline_target, sweep_indices)
        
        % Update parameter structure
        p_current = p;
        for i = 1:length(param_names)
            p_current.(param_names{i}) = x(i);
        end
        
        try
            % Calculonge Magic Formula output - FIXED: Use correct output variables
            [Fx_pred, ~, ~, ~, ~] = FXclass.Fx(p_current, longslip, slipangl, Fz, pressure, inclangl);
            
            % Calculonge base residuals
            residuals = Fx_pred - Fx_spline_target;

            % Normalize by sweep to ensure equal contribution
            unique_sweeps = unique(sweep_indices);
            for sweep_idx = unique_sweeps'
                sweep_mask = (sweep_indices == sweep_idx);
                if any(sweep_mask)
                    sweep_max_force = max(abs(Fx_spline_target(sweep_mask)));
                    if sweep_max_force > 0
                        residuals(sweep_mask) = residuals(sweep_mask) / sweep_max_force;
                    end
                end
            end
            
        catch ME
            fprintf('Error in Fx: %s\n', ME.message);
            residuals = 1000 * ones(size(Fx_spline_target));
        end
        end
    end

    methods

        function obj = optimisationParameters(obj)


            obj.params.params_to_optimize = {
    'PCX1', 'PDX1', 'PDX2', 'PDX3', ...  % Shape, Peak, and camber peak
    'PEX1', 'PEX2', 'PEX3', 'PEX4', ...  % Curvature factors  
    'PKX1', 'PKX2', 'PKX3', ...          % Stiffness factors
    'PHX1', 'PHX2', 'PVX1', 'PVX2', ...  % Shift factors
    'PPX1', 'PPX2', 'PPX3', 'PPX4', ...  % Pressure effects
};

            obj.params.options = optimoptions('lsqnonlin', ...
                    'Display', 'iter-detailed', ...
                    'MaxIterations', 500, ...
                    'MaxFunctionEvaluations', 200000, ...
                    'StepTolerance', 1e-8, ...
                    'FunctionTolerance', 1e-8, ...
                    'OptimalityTolerance', 1e-6, ...
                    'FiniteDifferenceType', 'central', ...
                    'Algorithm', 'levenberg-marquardt');
            
        end

        function obj = combinedoptimisationParameters(obj)
            obj.params.params_to_optimize = {
    'RBX1', 'RBX2', 'RBX3', ...     % Combined slip parameters
    'RCX1', 'REX1', 'REX2', ...     % Combined slip parameters
    'RHX1'                          % Combined slip parameters
};

            obj.params.options = optimoptions('lsqnonlin', ...
                    'Display', 'iter-detailed', ...
                    'MaxIterations', 2000, ...
                    'MaxFunctionEvaluations', 200000, ...
                    'StepTolerance', 1e-8, ...
                    'FunctionTolerance', 1e-8, ...
                    'OptimalityTolerance', 1e-6, ...
                    'FiniteDifferenceType', 'central', ...
                    'Algorithm', 'levenberg-marquardt');
        end

    end
end

