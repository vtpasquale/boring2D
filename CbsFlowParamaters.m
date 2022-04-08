classdef CbsFlowParamaters
    % Class for reading and storing CBS flow parameters (.par file extension)
    
    properties
        % Line 2
        restart % [logical]
        
        % Line 4
        energyCalculation  % [logical]
        
        % Line6
        Ux % [double] free stream value
        Uy % [double] free stream value
        P % [double] free stream value
        T % [double] free stream value
        
        % Line 8
        ntime % [integer] max number of pseudo time steps per real time step
        beta_opt % [logical] controls whether a constant beta (EPSILON) is applied throughout the domain, or locally varying
        epsilon % [double] small constant to avoid zero beta
        dtfixed % [double] fix local time step flag: 1 = ON, 0 = local, -1 = Global minimum of local values
        dtfix % [double] fixed local time step size
        iwrite % [double]
        
        % Line 10
        nRealTimesteps % [double]
        realTimestepSize % [double]
        
        % Line 12
        csafm % [double] time step safety factor
        theta1 % [double] time step parameter
        
        % Line 14
        Re % [double] Reynolds Number
        Pr % [double] Prandtl Number
        Ra % [double] Rayleigh Number
        Ri % [double] Richardson Number
        
        % Line 16
        convectionType % [double] (1 - Natural, 0 - Mixed/Forced)
        
        % Line 18 - steady state tolerances
        velocityFlag % [double] (0 - OFF, 1 - ON)
        velocityTol % [double] tolerance
        pressureFlag % [double] (0 - OFF, 1 - ON)
        pressureTol % [double] tolerance
        energyFlag % [double] (0 - OFF, 1 - ON)
        energyTol % [double] tolerance
        
        % Line 20 - output control
        paraviewFlag % [double] (0 - OFF, 1 - ON)
        tecplotFlag % [double] (0 - OFF, 1 - ON)
        localNusseltFlag % [double] (0 - OFF, 1 - ON)
        boundaryNusseltFlag % [double] ?
        
        % Line 22 - rumtime conrol updates
        runTimeControl % [double] (0 - OFF, 1 - ON)
    end
    
    methods
        function obj = CbsFlowParamaters(in)
            % Gmf class constructor
            if ischar(in)
                % check file extension
                [~,~,extension] = fileparts(in);
                switch extension
                    case '.par'
                        obj = obj.constructFromParFile(in);
                    otherwise
                        error('File extension %s not supported',extension)
                end
            else
                error('CbsFlowParamaters class constructor input type incorrect')
            end
        end % CbsFlowParamaters()
    end
    methods (Access=private)
        function obj = constructFromParFile(obj,in)
            % Construct from ascii file
            fid = fopen(in,'r');
                        
            skipLine = fgetl(fid);
            obj.restart = fscanf(fid,'%d\n',1); 
            
            skipLine = fgetl(fid);
            obj.energyCalculation = fscanf(fid,'%d\n',1); 
            
            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.Ux = str2num(splitLine{1}); 
            obj.Uy = str2num(splitLine{2}); 
            obj.P  = str2num(splitLine{3}); 
            obj.T  = str2num(splitLine{4}); 
            
            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.ntime = str2num(splitLine{1});
            obj.beta_opt = str2num(splitLine{2});
            obj.epsilon = str2num(splitLine{3});
            obj.dtfixed = str2num(splitLine{4});
            obj.dtfix = str2num(splitLine{5});
            obj.iwrite = str2num(splitLine{6});
            
            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.nRealTimesteps = str2num(splitLine{1});
            obj.realTimestepSize = str2num(splitLine{2});

            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.csafm = str2num(splitLine{1});
            obj.theta1 = str2num(splitLine{2});
            
            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.Re = str2num(splitLine{1});
            obj.Pr = str2num(splitLine{2});
            obj.Ra = str2num(splitLine{3});
            obj.Ri = str2num(splitLine{4});
            
            skipLine = fgetl(fid);
            obj.convectionType = fscanf(fid,'%d\n',1); 
            
            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.velocityFlag  = str2num(splitLine{1});
            obj.velocityTol   = str2num(splitLine{2});
            obj.pressureFlag  = str2num(splitLine{3});
            obj.pressureTol   = str2num(splitLine{4});
            obj.energyFlag    = str2num(splitLine{5});
            obj.energyTol     = str2num(splitLine{6});
            
            skipLine = fgetl(fid);
            getLine = fgetl(fid);
            splitLine = split(getLine);
            obj.paraviewFlag        = str2num(splitLine{1});
            obj.tecplotFlag         = str2num(splitLine{2});
            obj.localNusseltFlag    = str2num(splitLine{3});
            obj.boundaryNusseltFlag = str2num(splitLine{4});
            
            skipLine = fgetl(fid);
            obj.runTimeControl = fscanf(fid,'%d\n',1); 
            
            fclose(fid);
        end % constructFromPltFile()
    end
end

