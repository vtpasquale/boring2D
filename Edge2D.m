classdef Edge2D
    % Class for 2D boundary edges
    
    % Anthony Ricciardi
    % April 2022
    
    properties
        nodeIDs % [nEdge,2 uint32] node numbers
        x % [nEdge,2 double] x node locations
        y % [nEdge,2 double] y node locations
        boundaryID % [nEdge,1 uint32] boundary ID numbers
        
        gDof % [nEdge,4 uint32] global DOF from edge DOF. Assumes two dof per node in node order.
        
        n % [nEdge,2 double] edge unit normal vector
        length % [nEdge,1 double] edge length 
    end
    
    methods
        function obj = Edge2D(gmf)
            % Construct from Gamma Mesh object
            
            % Process node IDs
            [nEdge,mEdge] = size(gmf.edges);
            if nEdge < 1; error('No edges'); end
            if mEdge ~=3; error('mEdge~=3');     end
            if any(any(mod(gmf.edges,1)~=0)); error('Noninteger value'); end
            obj.nodeIDs = uint32(gmf.edges(:,1:2));
            obj.boundaryID = uint32(gmf.edges(:,3));
            
            % Process locations
            x = [ gmf.nodes(obj.nodeIDs(:,1),1),...
                  gmf.nodes(obj.nodeIDs(:,2),1)];
            y = [ gmf.nodes(obj.nodeIDs(:,1),2),...
                  gmf.nodes(obj.nodeIDs(:,2),2)];
            obj.x = x;
            obj.y = y;
            
            % global DOF index (assumes two dof per node in node order)
            obj.gDof =zeros(nEdge,4,'uint32');
            obj.gDof(:,[1,3]) = 2*obj.nodeIDs - 1;
            obj.gDof(:,[2,4]) = 2*obj.nodeIDs;
            
            % Length
            X12 = [x(:,2)-x(:,1),y(:,2)-y(:,1)];
            % obj.length = sqrt( (x(:,2)-x(:,1)).^2 + (y(:,2)-y(:,1)).^2 );
            obj.length = sqrt( sum(X12.^2,2) );
            obj.n = [X12(:,2), -X12(:,1)]./obj.length;
        end
    end
    
end

