classdef CbsRestartFile
    % Class for reading and storing CBS boundary definition (.var file extension)
    
    properties
        u1 % [nNodes,1 double]
        u2 % [nNodes,1 double]
        p  % [nNodes,1 double]
        T  % [nNodes,1 double]
    end
    
    methods
        function obj = CbsRestartFile(fileName,nNodes)
            % Gmf class constructor
            if ischar(fileName)
                % check file extension
                [~,~,extension] = fileparts(fileName);
                switch extension
                    case '.var'
                        obj = obj.constructFromVarFile(fileName,nNodes);
                    otherwise
                        error('File extension %s not supported',extension)
                end
            else
                error('CbsRestartFile class constructor input type incorrect')
            end
        end % CbsBoundaryDefinition()
    end
    methods (Access=private)
        function obj = constructFromVarFile(obj,in,nNodes)
            % Construct from ascii file
            fid = fopen(in,'r');
            
            skipLine = fgetl(fid);
            [data,count] = fscanf(fid,'%f',4*nNodes);
            if count~=4*nNodes; error('count~=4'); end
            obj.u1 = data(1:4:end);
            obj.u2 = data(2:4:end);
            obj.p  = data(3:4:end);
            obj.T  = data(4:4:end);
            
            fclose(fid);
        end % constructFromVarFile()
    end
end

