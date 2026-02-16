classdef FYclass < TyreModel & Parameters

    properties
        params % struct
    end

    methods(Static)
        function [Fy0,muy,dfz,Fz0_,dpi,By,Cy,Ey,Dy,Kya_,SVy,SHy] = Fy0(p,slipangl,Fz,pressure,inclangl)
            %MAGICFORMULA.V61.FY0
            % TODO: implement for turnslip (ZETAs):
            ZETA0 = 1;
            ZETA2 = 1;
            ZETA3 = 1;
            ZETA4 = 1;

            % (4.E4)
            gammaAst = sin(inclangl);
            gammaAst2 = gammaAst.^2;

            % (4.E1)
            Fz0_ = p.FNOMIN.*p.LFZO;

            % (4.E2a) 
            dfz = (Fz-Fz0_)./Fz0_;

            % (4.E2b) 
            pi0 = p.NOMPRES;
            dpi = (pressure-pi0)./pi0;
            dpi2 = dpi.^2;

            % (4.E21) 
            Cy = p.LCY.*p.PCY1;

            % (4.E23)
            % Note: Asterisk variant of LMUX is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            muy = (p.PDY1+p.PDY2.*dfz).*(1+p.PPY3.*dpi+p.PPY4.*dpi2).*(1-p.PDY3.*gammaAst2).*p.LMUY;

            % (4.E22) 
            Dy = muy.*Fz.*ZETA2;

            % (4.E25)
            Kya = p.PKY1.*Fz0_...
                .*(1+p.PPY1.*dpi)...
                .*(1-p.PKY3.*abs(gammaAst))...
                .*sin(p.PKY4.*atan(Fz./Fz0_./((p.PKY2+p.PKY5.*gammaAst2)...
                .*(1+p.PPY2.*dpi)))).*p.LKY.*ZETA3;

            % (4.E39); see also p.177 for explanation on eps/sign
            signKya = sign(Kya);
            signKya(~signKya) = 1;
            Kya_ = Kya + eps(Kya).*signKya;

            % (4.E28)
            % Note: Asterisk variant of LMUX is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            SVyg = ZETA2.*p.LKYC.*p.LMUY.*Fz.*(p.PVY3+p.PVY4.*dfz).*gammaAst;

            % (4.E30)
            Kyg0 = Fz.*(p.PKY6+p.PKY7.*dfz).*(1+p.PPY5.*dpi).*p.LKYC;

            % (4.E29)
            % Note: Asterisk variant of LMUX is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            SVy = ZETA2.*p.LMUY.*p.LVY.*Fz.*(p.PVY1+p.PVY2.*dfz)+SVyg;

            % (4.E27)
            SHy = p.LHY.*(p.PHY1+p.PHY2.*dfz)+(Kyg0.*gammaAst-SVyg)./Kya_*ZETA0+ZETA4-1;

            % (4.E20) 
            alphay = slipangl+SHy;
            alphaySgn = sign(slipangl);

            % (4.E24) 
            Ey = (p.PEY1+p.PEY2.*dfz)...
                .*(1+p.PEY5.*gammaAst.^2-(p.PEY3+p.PEY4.*gammaAst).*alphaySgn).*p.LEY;

            % (4.E26); see also p.177 for explanation on eps/sign
            signCy = sign(Cy);
            signCy(~signCy) = 1;
            By = Kya./(Cy.*Dy + eps(Cy).*signCy);

            % (4.E19)
            Fy0 = Dy.*sin(Cy.*atan(By.*alphay-Ey.*(By.*alphay-atan(By.*alphay))))+SVy;
        end

        function [Fy,muy,Gyk,Byk,Eyk,Fy0,gammaAst] = Fy(p,longslip,slipangl,Fz,pressure,inclangl)
        %MAGICFORMULA.V61.FY
        % TODO: implement for turnslip (ZETAs):
        ZETA2 = 1;
        
        [Fy0,muy,dfz,Fz0_,dpi,By,Cy,Ey,Dy,Kya_,SVy,SHy] = FYclass.Fy0(p,slipangl,Fz,pressure,inclangl);
        
        % (4.E4)
        gammaAst = sin(inclangl);
        
        % (4.E3)
        alphaAst = tan(slipangl);
        
        % (4.E62)
        Byk = p.LYKA.*(p.RBY1+p.RBY4*gammaAst.^2).*cos(atan(p.RBY2.*(alphaAst-p.RBY3)));
        
        % (4.E63)
        Cyk = p.RCY1;
        
        % (4.E67)
        DVyk = ZETA2*muy.*Fz.*(p.RVY1+p.RVY2.*dfz+p.RVY3*gammaAst).*cos(atan(p.RVY4.*alphaAst));
        
        % (4.E64)
        Eyk = p.REY1+p.REY2.*dfz;
        
        % (4.E65)
        SHyk = p.RHY1+p.RHY2.*dfz;
        
        % (4.E66)
        SVyk = p.LVYKA.*DVyk.*sin(p.RVY5.*atan(p.RVY6.*longslip));
        
        % (4.E61)
        kappaS = longslip + SHyk;
        
        % (4.E60)
        Gyk0 = cos(Cyk.*atan(Byk.*SHyk-Eyk.*(Byk.*SHyk-atan(Byk.*SHyk))));
        
        % (4.E59)
        Gyk = cos(Cyk.*atan(Byk.*kappaS-Eyk.*(Byk.*kappaS-atan(Byk.*kappaS))))./Gyk0;
        
        % (4.E58)
        Fy = Gyk.*Fy0 + SVyk;
        end


        %% Residual function with sweep normalization for Fy0
        function residuals = residual_fy0(x, p, param_names, ...
            slipangl, Fz, pressure, inclangl, Fy_spline_target, sweep_indices)
        
        % Update parameter structure
        p_current = p;
        for i = 1:length(param_names)
            p_current.(param_names{i}) = x(i);
        end
        
        try
            % Calculate Magic Formula output
            [Fy0_pred, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                FYclass.Fy0(p_current, slipangl, Fz, pressure, inclangl);
            
            % Calculate base residuals
            residuals = Fy0_pred - Fy_spline_target;
           
            % Normalize by sweep to ensure equal contribution
            unique_sweeps = unique(sweep_indices);
            for sweep_idx = unique_sweeps'
                sweep_mask = (sweep_indices == sweep_idx);
                if any(sweep_mask)
                    sweep_max_force = max(abs(Fy_spline_target(sweep_mask)));
                    if sweep_max_force > 0 
                        residuals(sweep_mask) = residuals(sweep_mask)/sweep_max_force; 
                
                    end
                end

            end
           
            
        catch ME
            fprintf('Error in Fy0: %s\n', ME.message);
            residuals = 1000 * ones(size(Fy_spline_target));
        end
        end

        %% Residual function with sweep normalization for combined
        function residuals = residual_fy(x, p, param_names, ...
            longslip, slipangl, Fz, pressure, inclangl, Fy_spline_target, sweep_indices)
        
        % Update parameter structure
        p_current = p;
        for i = 1:length(param_names)
            p_current.(param_names{i}) = x(i);
        end
        
        try
            % Calc Magic Formula output
            [Fy_pred, ~, ~, ~, ~, ~, ~] = FYclass.Fy(p_current, longslip, slipangl, Fz, pressure, inclangl);
            
            % Calculonge base residuals
            residuals = Fy_pred - Fy_spline_target;
           
            
            % Normalize by sweep to ensure equal contribution
            unique_sweeps = unique(sweep_indices);
            for sweep_idx = unique_sweeps'
                sweep_mask = (sweep_indices == sweep_idx);
                if any(sweep_mask)
                    sweep_max_force = max(abs(Fy_spline_target(sweep_mask)));
                    if sweep_max_force > 0
                        residuals(sweep_mask) = residuals(sweep_mask) / sweep_max_force;
                    end
                end
            end            
            
        catch ME
            fprintf('Error in Fy: %s\n', ME.message);
            residuals = 1000 * ones(size(Fy_spline_target));
        end
        end

end
    

    methods
            
        function obj = optimisationParameters(obj)
            obj.params = struct;
            obj.params.params_to_optimize = {
    'PCY1', 'PDY1', 'PDY2', 'PDY3', ...  % Shape, Peak, and camber peak
    'PEY1', 'PEY2', 'PEY3', 'PEY4', 'PEY5' ...  % Curvature factors  
    'PKY1', 'PKY2', 'PKY3', 'PKY4', 'PKY5' ...  % Stiffness factors
    'PKY6', 'PKY7', ...                  % Camber stiffness
    'PHY1', 'PHY2', ...                  % Horizontal shift
    'PVY1', 'PVY2', 'PVY3', 'PVY4', ...  % Vertical shift
    'PPY1', 'PPY2', 'PPY3', 'PPY4', 'PPY5' ... % Pressure effects
};
            obj.params.options = optimoptions('lsqnonlin', ...
                        'Display', 'iter-detailed', ...
                        'MaxIterations', 500, ...
                        'MaxFunctionEvaluations', 200000, ...
                        'StepTolerance', 1e-8, ...
                        'FunctionTolerance', 1e-8, ...
                        'OptimalityTolerance', 1e-6, ... 
                        'FiniteDifferenceType', 'central', ...
                        'Algorithm', 'interior-point');
        end

        function obj = combinedoptimisationParameters(obj)
            obj.params.params_to_optimize = {
    'RBY1', 'RBY2', 'RBY3', 'RBY4', ...  % Combined slip parameters
    'RCY1', 'REY1', 'REY2', ...          % Combined slip parameters
    'RVY1', 'RVY2', 'RVY3', 'RVY4', 'RVY5', 'RVY6', 'RHY1', 'RHY2'  % Combined slip parameters
};
            obj.params.options = optimoptions('lsqnonlin', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 200000, ...
    'StepTolerance', 1e-8, ...
    'FunctionTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-6, ...
    'FiniteDifferenceType', 'central', ...
    'Algorithm', 'interior-point');

        end 

    end
end
