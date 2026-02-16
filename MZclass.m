classdef MZclass < TyreModel & Parameters 

    properties
        params
    end

    methods(Static)
        function [Mz0,Bt,Ct,Dt,Et,Br,Cr,Dr,alphaCos,R0,alphar,alphat] = Mz0(p,slipangl,Fz,pressure,inclangl,Vcx)
            %MAGICFORMULA.V61.MZ0
            % TODO: implement for turnslip (ZETAs):
            ZETA0 = 1;
            ZETA2 = 1;
            ZETA5 = 1;
            ZETA6 = 1;
            ZETA7 = 1;
            ZETA8 = 1;
            
            % (4.E3)
            Vcy = -Vcx.*tan(slipangl);
            sgnVcx = sign(Vcx);
            alphaAst = tan(slipangl).*sgnVcx;
            
            % (4.E6a)
            Vc = sqrt(Vcx.^2 + Vcy.^2);
            Vc = Vc + eps(Vc);
            
            % (4.E6)
            alphaCos = Vcx./Vc;
            
            % (4.E4)
            gammaAst = sin(inclangl);
            gammaAst2 = gammaAst.^2;
            gammaAstAbs = abs(gammaAst);
            

            [Fy0,~,dfz,Fz0_,dpi,By,Cy,~,~,Kya_,SVy,SHy] = FYclass.Fy0(p,slipangl,Fz,pressure,0);
            dfz2 = dfz.^2;
            R0 = p.UNLOADED_RADIUS;
            
            % (4.E38), (4.E39)
            SHf = SHy + SVy./Kya_;
            
            % (4.E35)
            SHt = p.QHZ1 + p.QHZ2.*dfz + (p.QHZ3 + p.QHZ4.*dfz).*gammaAst;
            
            % (4.E34)
            alphat = alphaAst + SHt;
            
            % (4.E37)
            alphar = alphaAst + SHf;
            
            % (4.E42)
            Dt0 = Fz.*(R0./Fz0_).*(p.QDZ1+p.QDZ2.*dfz).*(1-p.PPZ1.*dpi).*p.LTR.*sgnVcx;
            
            % (4.E40)
            % Note: QBZ4 added although not included in Pacejka's 2012 book.
            % Note: QBZ6 eliminated because not included in 6.1.2 manual.
            % Note: Asterisk variant of LMUY is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            Bt = (p.QBZ1+p.QBZ2.*dfz+p.QBZ3.*dfz2)...
                .*(1+p.QBZ4.*gammaAst+p.QBZ5.*gammaAstAbs)...
                .*p.LKY./p.LMUY;
            
            % (4.E41)
            Ct = p.QCZ1;
            
            % (4.E43)
            Dt = Dt0.*(1+p.QDZ3.*gammaAstAbs+p.QDZ4.*gammaAst2).*ZETA5;
            
            % (4.E44)
            Et = (p.QEZ1+p.QEZ2.*dfz+p.QEZ3.*dfz2)...
                .*(1+(p.QEZ4+p.QEZ5.*gammaAst).*(2/pi).*atan(Bt.*Ct.*alphat));
            
            % (4.E45)
            % Note: Asterisk variant of LMUY is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            Br = (p.QBZ9.*p.LKY./p.LMUY+p.QBZ10.*By.*Cy).*ZETA6;
            
            % (4.E46)
            Cr = ZETA7;
            
            % (4.E47)
            % Note: Asterisk variant of LMUY is not used, because corresponding scaling
            %       factor is not defined in 6.1.2 manual.
            Dr = Fz.*R0.*((p.QDZ6+p.QDZ7.*dfz).*p.LRES.*ZETA2...
                +((p.QDZ8+p.QDZ9.*dfz).*(1+p.PPZ2.*dpi)...
                +(p.QDZ10+p.QDZ11.*dfz).*gammaAstAbs).*gammaAst.*p.LKZC.*ZETA0)...
                .*p.LMUY.*sgnVcx.*alphaCos + ZETA8 - 1;
            
            % (4.E33)
            t0 = Dt.*cos(Ct.*atan(Bt.*alphat-Et.*(Bt.*alphat-atan(Bt.*alphat)))).*alphaCos;
            
            % (4.E32)
            Mz0_ = -t0.*Fy0;
            
            % (4.E36)
            Mzr0 = Dr.*cos(Cr.*atan(Br.*alphar)).*alphaCos;
            
            % (4.E31)
            Mz0 = Mz0_ + Mzr0;
        end

        function Mz = Mz(p,longslip,slipangl,Fz,pressure,inclangl,Vcx, Fx, Fy)
            % Fx is final output for combined slip 
            % Fy is final output for combined slip 
            
            [~,Bt,Ct,Dt,Et,Br,Cr,Dr,alphaCos,R0,alphar,alphat] = MZclass.Mz0(p,slipangl,Fz,pressure,inclangl,Vcx);
            [~,~,~,~,~,~,Kxk] = FXclass.Fx0(p,longslip,Fz,pressure,inclangl);
            [~,~,dfz,Fz0_,~,~,~,~,~,Kya_] = FYclass.Fy0(p,slipangl,Fz,pressure,inclangl);
            [~, ~,Gyk,~,~,Fy0,gammaAst] = FYclass.Fy(p,longslip,slipangl,Fz,pressure,0);
            
            kappa2 = longslip.^2;
            KxkKya2 = (Kxk./Kya_).^2;
            
            % (4.E78)
            alphar2 = alphar.^2;
            alpharSgn = sign(alphar);
            alpharEq = sqrt(alphar2 + kappa2.*KxkKya2).*alpharSgn;
            
            % (4.E77)
            alphat2 = alphat.^2;
            alphatSgn = sign(alphat);
            alphatEq = sqrt(alphat2 + kappa2.*KxkKya2).*alphatSgn;
            % alphatEq = atan(sqrt(tan(alphat).^2+(Kxk./Kya_prime).^2.*kappa.^2)).*sign(alphat); % (A55)
            
            % (4.E76)
            s = R0.*(p.SSZ1+p.SSZ2.*(Fy/Fz0_)+(p.SSZ3+p.SSZ4.*dfz).*gammaAst).*p.LS;
            
            % (4.E75)
            Mzr = Dr.*cos(Cr.*atan(Br.*alpharEq)).*alphaCos;
            
            % (4.E74)
            Fy_ = Gyk.*Fy0;
            
            % (4.E73) Pneumatic trail
            t = Dt.*cos(Ct.*atan(Bt.*alphatEq-Et.*(Bt.*alphatEq-atan(Bt.*alphatEq)))).*alphaCos;
            
            % (4.E72)
            Mz_ = -t.*Fy_;
            
            % (4.E71)
            Mz = Mz_ + Mzr + s.*Fx;
        end



       %% Residual function with sweep normalization for MZ
        function residuals = residual_mz0(x, p, param_names, ...
        slipangl, Fz, pressure, inclangl, Vcx, Mz_spline_target, sweep_indices)
    
        % Update parameter structure
        p_current = p;
        for i = 1:length(param_names)
            p_current.(param_names{i}) = x(i);
        end
        
        try
            % Calculate Magic Formula output for MZ
            [Mz0_pred, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                    MZclass.Mz0(p_current, slipangl, Fz, pressure, inclangl, Vcx);
            
            % Calculate base residuals
            residuals = Mz0_pred - Mz_spline_target;
            
            % Normalize by sweep to ensure equal contribution
            unique_sweeps = unique(sweep_indices);
            for sweep_idx = unique_sweeps'
                sweep_mask = (sweep_indices == sweep_idx);
                if any(sweep_mask)
                    sweep_max_moment = max(abs(Mz_spline_target(sweep_mask)));
                    if sweep_max_moment > 0
                        residuals(sweep_mask) = residuals(sweep_mask) / sweep_max_moment;
                    end
                end
            end
            
        catch ME
            fprintf('Error in Mz0: %s\n', ME.message);
            residuals = 1000 * ones(size(Mz_spline_target));
        
        end
        
        end

        function residuals = residual_mz(x, p, param_names, longslip, ... 
        slipangl, Fz, pressure, inclangl, Vcx, Fx, Fy, Mz_spline_target, sweep_indices)

            % update parameter structure
            p_current = p;
            for i = 1:length(param_names)
                p_current.(param_names{i}) = x(i);
            end
            
            try
                % Calculate combined slip magic formula output for Mz
                Mz_pred = MZclass.Mz(p_current, longslip, slipangl, Fz, pressure, inclangl, Vcx, Fx, Fy);

                % Calculate the residuals
                residuals = Mz_pred - Mz_spline_target;

                % Normalise by sweep to ensure equal contribution 
                unique_sweeps = unique(sweep_indices);
                for sweep_idx = unique_sweeps'
                    sweep_mask = (sweep_indices == sweep_idx);
                    if any(sweep_mask) 
                        sweep_max_moment = max(abs(Mz_spline_target(sweep_mask))) ;
                        if sweep_max_moment > 0 
                            residuals(sweep_mask) = residuals(sweep_mask) / sweep_max_moment ;
                        end
                    end
                end
            catch ME
                fprintf(' Error in Mz: %s\n', ME.message);
                residuals = 1000* ones(size(Mz_spline_target));
            end
        end
    end

    methods

        function obj = optimisationParameters(obj)
            % Define which parameters to optimize for MZ
            obj.params.params_to_optimize = {
                'QBZ1', 'QBZ2', 'QBZ3', 'QBZ4', 'QBZ5', ...  % Pneumatic trail parameters
                'QBZ9', 'QBZ10', ...                         % Residual moment parameters
                'QCZ1', ...                                  % Pneumatic trail shape factor
                'QDZ1', 'QDZ2', 'QDZ3', 'QDZ4', 'QDZ6', 'QDZ7', 'QDZ8', 'QDZ9', 'QDZ10', 'QDZ11', ... % Peak factors
                'QEZ1', 'QEZ2', 'QEZ3', 'QEZ4', 'QEZ5', ...  % Curvature factors
                'QHZ1', 'QHZ2', 'QHZ3', 'QHZ4', ...          % Horizontal shift factors
                'PPZ1', 'PPZ2' ...                           % Pressure effects
            };
            % Optimization options
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
            % Define which parameters to optimise 
            obj.params.params_to_optimize = {
                'SSZ1', 'SSZ2', 'SSZ3', 'SSZ4',...             
                'LS'};                                  % Over turning moment scale factor 

            % Optimisation options
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

    end
end




