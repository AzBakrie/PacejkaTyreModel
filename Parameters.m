%% Parameters subclass
classdef (Abstract) Parameters
    %PARAMETERS Magic Formula 6.1.2 parameter set.
    % This class is a mixin: intended to be combined with another
    % superclass (like TireModel) via multiple inheritance.
    % 
    % Parameter names implemented according to:
    % https://functionbay.com/documentation/onlinehelp/Documents/Tire/MFTyre-MFSwift_Help.pdf
    %
    
    properties
        %% [MODEL]
        FITTYP
        TYRESIDE
        LONGVL
        VXLOW
        ROAD_INCREMENT
        ROAD_DIRECTION
        %% [DIMENSION]
        UNLOADED_RADIUS
        WIDTH
        RIM_RADIUS
        RIM_WIDTH
        ASPECT_RATIO
        %% [OPERATING_CONDITIONS]
        INFLPRES
        NOMPRES
        %% [INERTIA]
        MASS
        IXX
        IYY
        BELT_IXX
        BELT_IYY
        BELT_MASS
        GRAVITY
        %% [VERTICAL]
        FNOMIN
        VERTICAL_DAMPING
        VERTICAL_STIFFNESS
        MC_CONTOUR_A
        MC_CONTOUR_B
        BREFF
        DREFF
        FREFF
        Q_RE0
        Q_V1
        Q_V2
        Q_FZ2
        Q_FCX
        Q_FCY
        Q_FCY2
        Q_CAM
        Q_CAM1
        Q_CAM2
        Q_CAM3
        Q_FYS1
        Q_FYS2
        Q_FYS3
        PFZ1
        BOTTOM_OFFST
        BOTTOM_STIFF
        %% [STRUCTURAL]
        LONGITUDINAL_STIFFNESS
        LATERAL_STIFFNESS
        YAW_STIFFNESS
        FREQ_LAT
        FREQ_LONG
        FREQ_WINDUP
        FREQ_YAW
        DAMP_LAT
        DAMP_LONG
        DAMP_RESIDUAL
        DAMP_VLOW
        DAMP_WINDUP
        DAMP_YAW
        Q_BVX
        Q_BVT
        PCFX1
        PCFX2
        PCFX3
        PCFY1
        PCFY2
        PCFY3
        PCMZ1
        %% [CONTACT_PATCH]
        Q_RA1
        Q_RA2
        Q_RB1
        Q_RB2
        ELLIPS_SHIFT
        ELLIPS_LENGTH
        ELLIPS_HEIGHT
        ELLIPS_ORDER
        ELLIPS_MAX_STEP
        ELLIPS_NWIDTH
        ELLIPS_NLENGTH
        ENV_C1
        ENV_C2
        Q_A1
        Q_A2
        %% [INFLATION_PRESSURE_RANGE]
        PRESMIN
        PRESMAX
        %% [VERTICAL_FORCE_RANGE]
        FZMAX
        FZMIN
        %% [LONG_SLIP_RANGE]
        KPUMIN
        KPUMAX
        %% [SLIP_ANGLE_RANGE]
        ALPMIN
        ALPMAX
        %% [INCLINATION_ANGLE_RANGE]
        CAMMIN
        CAMMAX
        %% [SCALING_COEFFICIENTS]
        LFZO
        LCX
        LMUX
        LEX
        LKX
        LHX
        LVX
        LCY
        LMUY
        LEY
        LKY
        LKYC
        LKZC
        LHY
        LVY
        LTR
        LRES
        LXAL
        LYKA
        LVYKA
        LS
        LMX
        LVMX
        LMY
        LMP
        %% [LONGITUDINAL_COEFFICIENTS]
        PCX1
        PDX1
        PDX2
        PDX3
        PEX1
        PEX2
        PEX3
        PEX4
        PKX1
        PKX2
        PKX3
        PHX1
        PHX2
        PVX1
        PVX2
        RBX1
        RBX2
        RBX3
        RCX1
        REX1
        REX2
        RHX1
        PPX1
        PPX2
        PPX3
        PPX4
        %% [OVERTURNING_COEFFICIENTS]
        QSX1
        QSX2
        QSX3
        QSX4
        QSX5
        QSX6
        QSX7
        QSX8
        QSX9
        QSX10
        QSX11
        QSX12
        QSX13
        QSX14
        PPMX1
        %% [LATERAL_COEFFICIENTS]
        PCY1
        PDY1
        PDY2
        PDY3
        PEY1
        PEY2
        PEY3
        PEY4
        PEY5
        PKY1
        PKY2
        PKY3
        PKY4
        PKY5
        PKY6
        PKY7
        PHY1
        PHY2
        PVY1
        PVY2
        PVY3
        PVY4
        RBY1
        RBY2
        RBY3
        RBY4
        RCY1
        REY1
        REY2
        RHY1
        RHY2
        RVY1
        RVY2
        RVY3
        RVY4
        RVY5
        RVY6
        PPY1
        PPY2
        PPY3
        PPY4
        PPY5
        %% [ROLLING_COEFFICIENTS]
        QSY1
        QSY2
        QSY3
        QSY4
        QSY5
        QSY6
        QSY7
        QSY8
        %% [ALIGNING_COEFFICIENTS]
        QBZ1
        QBZ2
        QBZ3
        QBZ4
        QBZ5
        QBZ9
        QBZ10
        QCZ1
        QDZ1
        QDZ2
        QDZ3
        QDZ4
        QDZ6
        QDZ7
        QDZ8
        QDZ9
        QDZ10
        QDZ11
        QEZ1
        QEZ2
        QEZ3
        QEZ4
        QEZ5
        QHZ1
        QHZ2
        QHZ3
        QHZ4
        SSZ1
        SSZ2
        SSZ3
        SSZ4
        PPZ1
        PPZ2
        %% [TURNSLIP_COEFFICIENTS]
        PDXP1
        PDXP2
        PDXP3
        PKYP1
        PDYP1
        PDYP2
        PDYP3
        PDYP4
        PHYP1
        PHYP2
        PHYP3
        PHYP4
        PECP1
        PECP2
        QDTP1
        QCRP1
        QCRP2
        QBRP1
        QDRP1
    end

    methods
        function obj = initialGuess(obj)

            % Initialize
            obj.params.x0 = [];
            obj.params.lb = [];
            obj.params.ub = [];
            % Nominal values
            obj.p.FNOMIN = 150 * 4.44822;
            obj.p.NOMPRES = 97000;
            % Scaling factors
            obj.p.LFZO = 1.0; obj.p.LCY = 1.0; obj.p.LMUY = 1.0; obj.p.LKY = 1.0;
            obj.p.LHY = 1.0; obj.p.LVY = 1.0; obj.p.LKYG = 1.0; obj.p.LKYC = 1.0; obj.p.LEY = 1.0;

            %% FY VALUES


            % Main parameters - initial guesses
            obj.p.PCY1 = 1.5;   obj.p.PDY1 = 1.7;   obj.p.PDY2 = -0.1;   obj.p.PDY3 = 20.0;
            obj.p.PEY1 = 1.5;   obj.p.PEY2 = 0.0;   obj.p.PEY3 = 0.0;    obj.p.PEY4 = 0.0;  obj.p.PEY5 = 0.0;
            obj.p.PKY1 = 2.5;   obj.p.PKY2 = 0.0;   obj.p.PKY3 = -1.0;   obj.p.PKY4 = 0.3;  obj.p.PKY5 = 1.0;
            obj.p.PKY6 = 3.0;   obj.p.PKY7 = 4.0;   obj.p.PHY1 = 0.0;    obj.p.PHY2 = 0.0;
            obj.p.PVY1 = 0.0;   obj.p.PVY2 = 0.0;   obj.p.PVY3 = 0.0;    obj.p.PVY4 = 0.0;
            obj.p.PPY1 = -10.0; obj.p.PPY2 = 10.0;  obj.p.PPY3 = 3.0;    obj.p.PPY4 = 0.0; obj.p.PPY5 = 0.0;
    
            % Scaling factors
            obj.p.LYKA = 1.0; obj.p.LVYKA = 1.0; obj.p.LHY = 1.0;
            
            % Main parameters - initial guesses (combined) 
            obj.p.RBY1 = 10.0;    obj.p.RBY2 = 8.0;    obj.p.RBY3 = 0.001;  obj.p.RBY4 = 1.0;
            obj.p.RCY1 = 1.2;     obj.p.REY1 = -0.5;   obj.p.REY2 = 0.2;
            obj.p.RVY1 = 0.01;    obj.p.RVY2 = 0.005;  obj.p.RVY3 = 0.1;    obj.p.RVY4 = 10.0;
            obj.p.RVY5 = 1.5;     obj.p.RVY6 = 0.8;    obj.p.RHY1 = 0.001;  obj.p.RHY2 = 0.0005;

        %% FX VALUES

            % Scaling factors
            obj.p.LCX = 1.0; obj.p.LMUX = 1.0; obj.p.LKX = 1.0;
            obj.p.LHX = 1.0; obj.p.LVX = 1.0; obj.p.LEX = 1.0;
    
            % Main parameters - initial guesses
            obj.p.PCX1 = 1.5;   obj.p.PDX1 = 2.0;   obj.p.PDX2 = -0.1;   obj.p.PDX3 = 20.0;
            obj.p.PEX1 = 1.5;   obj.p.PEX2 = 0.0;   obj.p.PEX3 = 0.0;    obj.p.PEX4 = 0.0;
            obj.p.PKX1 = 2.5;   obj.p.PKX2 = 0.0;   obj.p.PKX3 = -1.0;
            obj.p.PHX1 = 0.0;   obj.p.PHX2 = 0.0;   obj.p.PVX1 = 0.0;    obj.p.PVX2 = 0.0;
            obj.p.PPX1 = -10.0; obj.p.PPX2 = 10.0;  obj.p.PPX3 = 3.0;    obj.p.PPX4 = 0.0;

            % Scaling factors
            obj.p.LXAL = 1.0;
            
            % Main parameters - initial guesses
            obj.p.RBX1 = 10.0;    obj.p.RBX2 = 8.0;    obj.p.RBX3 = 0.001;
            obj.p.RCX1 = 1.2;     obj.p.REX1 = -0.5;   obj.p.REX2 = 0.2;
            obj.p.RHX1 = 0.001;



        %% MZ VALUES
        % Nominal values
            obj.p.UNLOADED_RADIUS = 0.3;  % Typical unloaded radius in meters
            
            % Scaling factors for MZ
            obj.p.LTR = 1.0; obj.p.LRES = 1.0; obj.p.LKZC = 1.0;    obj.p.LS = 1.0;
            
            % Main parameters for MZ - initial guesses
            obj.p.QBZ1 = 10.0;    obj.p.QBZ2 = 2.0;     obj.p.QBZ3 = 1.0;     obj.p.QBZ4 = 0.0;     obj.p.QBZ5 = -0.5;
            obj.p.QBZ9 = 10.0;    obj.p.QBZ10 = 0.0;
            obj.p.QCZ1 = 1.5;
            obj.p.QDZ1 = 0.1;     obj.p.QDZ2 = 0.0;     obj.p.QDZ3 = 0.0;     obj.p.QDZ4 = 0.0;     obj.p.QDZ6 = 0.0;
            obj.p.QDZ7 = 0.0;     obj.p.QDZ8 = 0.0;     obj.p.QDZ9 = 0.0;     obj.p.QDZ10 = 0.0;    obj.p.QDZ11 = 0.0;
            obj.p.QEZ1 = 0.5;     obj.p.QEZ2 = 0.0;     obj.p.QEZ3 = 0.0;     obj.p.QEZ4 = 0.0;     obj.p.QEZ5 = 0.0;
            obj.p.QHZ1 = 0.0;     obj.p.QHZ2 = 0.0;     obj.p.QHZ3 = 0.0;     obj.p.QHZ4 = 0.0;
            obj.p.PPZ1 = 0.0;     obj.p.PPZ2 = 0.0;

       % combined slip MZ parameters
            obj.p.SSZ1 = 1.0;    obj.p.SSZ2= 0.5;   obj.p.SSZ3 = 0.0;   obj.p.SSZ4 = 0.0;      


            %% FY BOUNDS
            % Define bounds in a struct
            if isa(obj, 'FYclass')
                bounds = struct( ...
                    'PCY1',  [1.0, 3.0], ...
                    'PDY1',  [0.8, 4.0], ...
                    'PDY2',  [-1.0, 50.0], ...
                    'PDY3',  [-1.0, 50.0], ...
                    'PEY1',  [1.0, 2.5], ...
                    'PEY2',  [-10.0, 10.0], ...
                    'PEY3',  [-10.0, 10.0], ...
                    'PEY4',  [-10.0, 10.0], ...
                    'PEY5',  [-1.0, 1.0], ...
                    'PKY1',  [0.0, 50.0], ...
                    'PKY2',  [0.0, 2.0], ...
                    'PKY3',  [-20.0, 20.0], ...
                    'PKY4',  [0.0, 2.0], ...
                    'PKY5',  [0.0, 5.0], ...
                    'PKY6',  [-20.0, 20.0], ...
                    'PKY7',  [-20.0, 20.0], ...
                    'PHY1',  [-0.2, 0.2], ...
                    'PHY2',  [-0.2, 0.2], ...
                    'PVY1',  [-0.2, 0.2], ...
                    'PVY2',  [-0.2, 0.2], ...
                    'PVY3',  [-10.0, 10.0], ...
                    'PVY4',  [-10.0, 10.0], ...
                    'PPY1',  [-500.0, 500.0], ...
                    'PPY2',  [-500.0, 500.0], ...
                    'PPY3',  [-500.0, 500.0], ...
                    'PPY4',  [-500.0, 500.0], ...
                    'PPY5',  [-500.0, 500.0], ...
                    'RBY1',  [-50.0, 100.0], ...
                    'RBY2',  [-50.0, 100.0], ...
                    'RBY3',  [-50.0, 100.0], ...
                    'RBY4',  [0.0, 200.0], ...
                    'RCY1',  [-10.0, 10.0], ...
                    'REY1',  [-10.0, 10.0], ...
                    'REY2',  [-10.0, 10.0], ...
                    'RVY1',  [-30.0, 30.0], ...
                    'RVY2',  [-30.0, 30.0], ...
                    'RVY3',  [-30.0, 30.0], ...
                    'RVY4',  [-30.0, 30.0], ...
                    'RVY5',  [-30.0, 30.0], ...
                    'RVY6',  [-30.0, 30.0], ...
                    'RHY1',  [-30.0, 30.0], ...
                    'RHY2',  [-30.0, 30.0] ...
                );
            elseif isa(obj, 'FXclass')
                bounds = struct( ...
                'PCX1', [0.8, 4.0], ...
                'PDX1', [1.0, 4.0], ...
                'PDX2', [-1.0, 50.0], ...
                'PDX3', [-1.0, 50.0], ...
                'PEX1', [0.0, 2.5], ...
                'PEX2', [-10.0, 10.0], ...
                'PEX3', [-10.0, 10.0], ...
                'PEX4', [-10.0, 10.0], ...
                'PKX1', [0.0, 100.0], ...
                'PKX2', [0.0, 2.0], ...
                'PKX3', [-20.0, 20.0], ...
                'PHX1', [-0.2, 0.2], ...
                'PHX2', [-0.2, 0.2], ...
                'PVX1', [-0.2, 0.2], ...
                'PVX2', [-0.2, 0.2], ...
                'PPX1', [-500.0, 500.0], ...
                'PPX2', [-500.0, 500.0], ...
                'PPX3', [-500.0, 500.0], ...
                'PPX4', [-500.0, 500.0], ...
                'RBX1', [0.0, 150.0], ...
                'RBX2', [0.0, 150.0], ...
                'RBX3', [-0.5, 0.5], ...
                'RCX1', [0.5, 2.0], ...
                'REX1', [-2.0, 0.0], ...
                'REX2', [-1.0, 1.0], ...
                'RHX1', [-0.1, 0.1] ...
            );
            elseif isa(obj, 'MZclass')
                bounds = struct( ...
                'QBZ1', [0.0, 50.0], ...
                'QBZ9', [0.0, 50.0], ...
                'QBZ2', [-10.0, 10.0], ...
                'QBZ3', [-10.0, 10.0], ...
                'QBZ4', [-10.0, 10.0], ...
                'QBZ5', [-10.0, 10.0], ...
                'QBZ10', [-10.0, 10.0], ...
                'QCZ1', [0.5, 2.5], ...
                'QDZ1', [-1.0, 1.0], ...
                'QDZ2', [-1.0, 1.0], ...
                'QDZ3', [-1.0, 1.0], ...
                'QDZ4', [-1.0, 1.0], ...
                'QDZ6', [-1.0, 1.0], ...
                'QDZ7', [-1.0, 1.0], ...
                'QDZ8', [-1.0, 1.0], ...
                'QDZ9', [-1.0, 1.0], ...
                'QDZ10', [-1.0, 1.0], ...
                'QDZ11', [-1.0, 1.0], ...
                'QEZ1', [-5.0, 5.0], ...
                'QEZ2', [-5.0, 5.0], ...
                'QEZ3', [-5.0, 5.0], ...
                'QEZ4', [-5.0, 5.0], ...
                'QEZ5', [-5.0, 5.0], ...
                'QHZ1', [-0.1, 0.1], ...
                'QHZ2', [-0.1, 0.1], ...
                'QHZ3', [-0.1, 0.1], ...
                'QHZ4', [-0.1, 0.1], ...
                'PPZ1', [-10.0, 10.0], ...
                'PPZ2', [-10.0, 10.0], ...
                'SSZ1', [-10.0, 3.0], ...
                'SSZ2', [-2.0, 2.0], ...
                'SSZ3', [-1.0, 1.0], ...
                'SSZ4', [-1.0, 1.0], ...
                'LS', [0.0, 3.0] ...
            );
            end


           % Loop over parameters
            for i = 1:length(obj.params.params_to_optimize)
                param_name = obj.params.params_to_optimize{i};
        
                % Current value
                obj.params.x0 = [obj.params.x0; obj.p.(param_name)];
        
                % Bounds
                b = bounds.(param_name);
                obj.params.lb = [obj.params.lb; b(1)];
                obj.params.ub = [obj.params.ub; b(2)];
            end 

        end


    end
end
