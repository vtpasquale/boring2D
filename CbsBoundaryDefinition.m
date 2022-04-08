classdef CbsBoundaryDefinition
    % Class for reading and storing CBS boundary definition (.bco file extension)
    
%     Boundary Condition File
%     -------------------------
%     The 1st line defines how many boundary side types are present in the mesh (NFLAG). Following this will then be NFLAG
%     entries in the FLAG_LIST array. This array replaces boundary side numbers in the PLT file with boundary flag codes.
%     This procedure is to avoid regenerating the mesh each time you change the boundary conditions. Currently, The available
%     flag codes are:
%     
%     !     500 - adiabatic with prescribed velocity
%     !     501 - constant temperature with no-slip (T = 1)
%     !     502 - constant temperature with no-slip(T = 0)
%     !     503 - constant temperature (T = 0) prescribed velocity (u=1,v=0)
%     !     504 - pressure boundary
%     !     506 - velocity symmetry, no-flux energy
%     !     507 - Backward Facing Step (Re=229) Parabolic Boundary (with v=0)
    
    properties
        flagList % [nFlags,1 double] corresponds to boundary side numbers in the PLT file
        flagCode % [nFlags,1 double] boundary codes
    end
    
    methods
        function obj = CbsBoundaryDefinition(in)
            % Gmf class constructor
            if ischar(in)
                % check file extension
                [~,~,extension] = fileparts(in);
                switch extension
                    case '.bco'
                        obj = obj.constructFromBcoFile(in);
                    otherwise
                        error('File extension %s not supported',extension)
                end
            else
                error('CbsBoundaryDefinition class constructor input type incorrect')
            end
        end % CbsBoundaryDefinition()
    end
    methods (Access=private)
        function obj = constructFromBcoFile(obj,in)
            % Construct from ascii file
            fid = fopen(in,'r');
            nFlag = fscanf(fid,'%d\n',1); 
            [flags,count] = fscanf(fid,'%d\n',2*nFlag);
            if count~=2*nFlag; error('count~=2*nFlag'); end
            obj.flagList = flags(1:2:end);
            obj.flagCode = flags(2:2:end);
            fclose(fid);
        end % constructFromBcoFile()
    end
end

