classdef Gmf
    % Class to store and process Gamma Mesh Format developed at INRIA.
    % 
    % Useful information:
    % https://github.com/LoicMarechal/libMeshb/
    % https://pyamg.saclay.inria.fr/download/vizir/vizir4_user_guide.pdf
    %
    % This class is intended only to provide processing for ascii format
    % files and simple 2D meshes. 
    
    % Anthony Ricciardi

    properties      
        nodes % [nNodes,2 double] Node locations (X,Y).
        edges % [nNodes,3 double] Edge nodes numbers and boundary number. 
        tri   % [nSurfTrias,3 double] Node indices for triangular boundary surface faces.        
    end   
    methods
        function obj = Gmf(in)
            % Gmf class constructor
            if ischar(in)
                % check file extension
                [~,~,extension] = fileparts(in);
                switch extension
                    case '.mesh'
                        obj = obj.constructFromMeshFile(in);
                    case '.plt'
                        obj = obj.constructFromPltFile(in);
                    otherwise
                        error('Mesh file extension %s not supported',extension)
                end
            else
                error('Gmf class constructor input type incorrect')
            end
        end % Gmf()
        
        function writeVTK(obj,filename)
            % Write the mesh to a .vtk file.
            % https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
            fid = fopen(filename,'w+');
            fprintf(fid,'# vtk DataFile Version 3.0\n');
            fprintf(fid,'My example\n');
            fprintf(fid,'ASCII\n');
            
            fprintf(fid,'\nDATASET UNSTRUCTURED_GRID\n');
            fprintf(fid,'POINTS %d float\n',size(obj.nodes,1));
            fprintf(fid,'%f %f 0\n',obj.nodes');
            
            nElements = size(obj.tri,1);
            nCellData = 4*nElements;
            fprintf(fid,'CELLS %d %d\n',nElements,nCellData);
            fprintf(fid,'3 %d %d %d\n',obj.tri'-1);
%             if obj.nSurfTrias; fprintf(fid,'3 %d %d %d\n',               obj.tri'-1);   end
%             if obj.nSurfQuads; fprintf(fid,'4 %d %d %d %d\n',            obj.quad'-1);  end
%             if obj.nVolTets;   fprintf(fid,'4 %d %d %d %d\n',            obj.tet'-1);   end
%             if obj.nVolPents5; fprintf(fid,'5 %d %d %d %d %d\n',         obj.pent5'-1); end
%             if obj.nVolPents6; fprintf(fid,'6 %d %d %d %d %d %d\n',      obj.pent6'-1); end
%             if obj.nVolHexs;   fprintf(fid,'8 %d %d %d %d %d %d %d %d\n',obj.hex'-1);   end
            
            fprintf(fid,'CELL_TYPES %d\n',nElements);
            fprintf(fid,'%d\n',5*ones(nElements,1));
%             fprintf(fid,'%d\n',[...
%                 5* ones(obj.nSurfTrias, 1);
%                 9* ones(obj.nSurfQuads, 1);
%                 10*ones(obj.nVolTets,   1);
%                 14*ones(obj.nVolPents5, 1);
%                 13*ones(obj.nVolPents6, 1);
%                 12*ones(obj.nVolHexs,   1)  ]);
            
            fclose(fid);
        end % writeVTK()
    end 
    
    
    methods (Access=private)
        function obj = constructFromMeshFile(obj,in)
            % Construct from Gamma Mesh Format ascii file
            fid = fopen(in,'r');
                        
            % Process header
            meshFormatText = fscanf(fid,'%s',1);
            if ~strcmpi(meshFormatText,'MeshVersionFormatted')
                error('File header format issue')
            end
            meshVersionFormat = fscanf(fid,'%d',1);
            if meshVersionFormat~=3
                error('MeshVersionFormatted~=3')
            end
            dimensionText = fscanf(fid,'%s',1);
            if ~strcmpi(dimensionText,'Dimension')
                error('File header format issue')
            end
            dimensionSize = fscanf(fid,'%d',1);
            if dimensionSize~=2
                error('Dimension~=2')
            end
            
            % Read nodes
            nodesText = fscanf(fid,'%s',1);
            if ~strcmpi(nodesText,'Vertices')
                error('Issue locating Vertices data block.')
            end
            nNodes = fscanf(fid,'%d',1);
            [nodes0,count] = fscanf(fid,'%f',3*nNodes);
            if count~=(3*nNodes); error('Issue processing Vertices'); end
            if any( nodes0(3:3:end) ~= 0 ); error('Some vertices have nonzero Z value.'); end
            obj.nodes = [nodes0(1:3:end),nodes0(2:3:end)];
            clear nodes0
            
            % Read corners
            cornersText = fscanf(fid,'%s',1);
            if ~strcmpi(cornersText,'Corners')
                error('Issue locating Corners data block.')
            end
            nCorners = fscanf(fid,'%d',1);
            if nCorners~=0; error('nCorners~=0'); end
            
            % Read edges
            edgesText = fscanf(fid,'%s',1);
            if ~strcmpi(edgesText,'Edges')
                error('Issue locating Edges data block.')
            end
            nEdges = fscanf(fid,'%d',1);            
            [edges0,count] = fscanf(fid,'%d',3*nEdges);
            if count~=(3*edges0); error('Issue processing Edges'); end
            obj.edges = [edges0(1:3:end),edges0(2:3:end),edges0(3:3:end)];
            clear edges0
            
            % Read EdgesP1Ordering
            edgesP1OrderingText = fscanf(fid,'%s',1);
            if ~strcmpi(edgesP1OrderingText,'EdgesP1Ordering')
                error('Issue locating EdgesP1Ordering data block.')
            end
            fscanf(fid,'%d',3); % skip data
            
            % Read triangles
            trianglesText = fscanf(fid,'%s',1);
            if ~strcmpi(trianglesText,'Triangles')
                error('Issue locating Triangles data block.')
            end
            nTriangles = fscanf(fid,'%d',1); 
            [tri0,count] = fscanf(fid,'%d',4*nTriangles);
            if count ~=(4*nTriangles); error('Processing issue'); end
            obj.tri = [tri0(1:4:end),tri0(2:4:end),tri0(3:4:end)];
            clear tri0
            
            % Read TrianglesP1Ordering
            trianglesP1OrderingText = fscanf(fid,'%s',1);
            if ~strcmpi(trianglesP1OrderingText,'TrianglesP1Ordering')
                error('Issue locating TrianglesP1Ordering data block.')
            end
            fscanf(fid,'%d',10); % skip data
            
            % Confirm end text
            endText = fscanf(fid,'%s',1);
            if ~strcmpi(endText,'end')
                error('Issue locating End marker.')
            end
            
            fclose(fid);
        end % constructFromMeshFile()
        
        function obj = constructFromPltFile(obj,in)
            % Construct from CBSFlow ascii file
            fid = fopen(in,'r');
                        
            % Process header
            nTriangles = fscanf(fid,'%d',1);
            nNodes = fscanf(fid,'%d',1);
            nEdges = fscanf(fid,'%d',1);

            % Read triangles
            [tri0,count] = fscanf(fid,'%d',4*nTriangles);
            if count ~=(4*nTriangles); error('Processing issue'); end
            obj.tri = [tri0(2:4:end),tri0(3:4:end),tri0(4:4:end)];
            elementNumber = tri0(1:4:end);
            if any(elementNumber~=(1:nTriangles).')
                error('Elements not in order - add sorting.')
            end
            clear tri0 elementNumber
            
            % Read nodes
            [nodes0,count] = fscanf(fid,'%f',3*nNodes);
            if count~=(3*nNodes); error('Issue processing Nodes'); end
            obj.nodes = [nodes0(2:3:end),nodes0(3:3:end)];
            nodeNumber = nodes0(1:3:end);
            if any(nodeNumber~=(1:nNodes).')
                error('Elements not in order - add sorting.')
            end
            clear nodes0 nodeNumber
            
            % Read edges            
            [edges0,count] = fscanf(fid,'%d',4*nEdges);
            if count~=(4*edges0); error('Issue processing Edges'); end
            obj.edges = [edges0(1:4:end),edges0(2:4:end),edges0(4:4:end)];
            clear edges0
            
            fclose(fid);
        end % constructFromPltFile()
        
        
    end 
    
end
