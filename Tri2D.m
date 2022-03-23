classdef Tri2D
    % Class for 2D triangular elements
    
    % Anthony Ricciardi
    % March 2022
    
    properties
        nodeIDs % [3,nTri uint32] node numbers
        x % [3,nTri double] x node locations
        y % [3,nTri double] y node locations
        

        invJ % [2,2,nTri Matrix3D] inverse of Jacobian matrix (constant inside elements)
        area % [:,:,nTri double] element area
        dNdX % [2,3,nTri Matrix3D] physical shape function derivatives (constant inside elements)
        B    % [3,6,nTri Matrix3D] velocity-strain rate matrix (constant inside elements)
        
        M    % [6,6,nTri double] mass matrix
        K    % [6,6,nTri double] stiffness matrix
    end
    
    methods
        function obj = Tri2D(gmf)
            % Construct from Gamma Mesh object
            
            % Process node IDs
            [nTri,mTri] = size(gmf.tri);
            if nTri < 1; error('No elements'); end
            if mTri ~=3; error('mTri~=3');     end
            if any(any(mod(gmf.tri,1)~=0)); error('Noninteger node number'); end
            obj.nodeIDs = uint32(gmf.tri);
            
            % Process locations
            x = [ gmf.nodes(obj.nodeIDs(:,1),1),...
                  gmf.nodes(obj.nodeIDs(:,2),1),...
                  gmf.nodes(obj.nodeIDs(:,3),1)];
            y = [ gmf.nodes(obj.nodeIDs(:,1),2),...
                  gmf.nodes(obj.nodeIDs(:,2),2),...
                  gmf.nodes(obj.nodeIDs(:,3),2)];
            obj.x = x;
            obj.y = y;
            
            % Process Jacobian - constant inside each element
            detJ = (x(:,2)-x(:,1)).*(y(:,3)-y(:,1)) - (x(:,3)-x(:,1)).*(y(:,2)-y(:,1));
            invJ = zeros(2,2,nTri);
            invDet = zeros(1,1,nTri);
            invDet(1,1,:) = (1./detJ);
            invJ(1,1,:) = y(:,3)-y(:,1);
            invJ(1,2,:) = -1*(y(:,2)-y(:,1));
            invJ(2,1,:) = -1*(x(:,3)-x(:,1));
            invJ(2,2,:) = (x(:,2)-x(:,1));
            invJ = invDet.*invJ;
            invJ = Matrix3D(invJ);
            obj.invJ = invJ;
            obj.area = zeros(1,1,nTri);
            obj.area(1,1,:) = 0.5*detJ;
            
            % Velocity-strain rate matrix - constant inside each element
            % CMPW (7.2-8) fast analytic method (membrane strain-dispacement)
            B = zeros(3,6,nTri);
            B(1,1,:) = y(:,2)-y(:,3);
            B(1,3,:) = y(:,3)-y(:,1);
            B(1,5,:) = y(:,1)-y(:,2);
            B(2,2,:) = x(:,3)-x(:,2);
            B(2,4,:) = x(:,1)-x(:,3);
            B(2,6,:) = x(:,2)-x(:,1);
            B(3,1,:) = x(:,3)-x(:,2);
            B(3,2,:) = y(:,2)-y(:,3);
            B(3,3,:) = x(:,1)-x(:,3);
            B(3,4,:) = y(:,3)-y(:,1);
            B(3,5,:) = x(:,2)-x(:,1);
            B(3,6,:) = y(:,1)-y(:,2);
            B = invDet.*B;
            obj.B = Matrix3D(B);
            
            % Shape function derivatives are constant inside the element
            dNdxi  = [-1, 1, 0];
            dNdeta = [-1, 0, 1];
            dNXi = [dNdxi;dNdeta];
            dNXiNd = repmat(dNXi,[1,1,nTri]);
            
            % Physical shape function derivatives
            obj.dNdX = invJ*dNXiNd;
            
%             % Alternate B computation
%             B2 = zeros(3,6,nTri);
%             B2(1,1,:) = obj.dNdX(1,1,:);
%             B2(1,3,:) = obj.dNdX(1,2,:);
%             B2(1,5,:) = obj.dNdX(1,3,:);
%             B2(2,2,:) = obj.dNdX(2,1,:);
%             B2(2,4,:) = obj.dNdX(2,2,:);
%             B2(2,6,:) = obj.dNdX(2,3,:);
%             B2(3,1,:) = obj.dNdX(2,1,:);
%             B2(3,2,:) = obj.dNdX(1,1,:);
%             B2(3,3,:) = obj.dNdX(2,2,:);
%             B2(3,4,:) = obj.dNdX(1,2,:);
%             B2(3,5,:) = obj.dNdX(2,3,:);
%             B2(3,6,:) = obj.dNdX(1,3,:);
%             B-B2 % check == small
            
        end
        
        function obj = computeMass(obj)
            nTri = size(obj.nodeIDs,1);
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            % shapeFunctions = @(r,s) [1-r-s, r, s];
            N1 = [1-r(1)-s(1),           0, r(1),    0, s(1),   0   ;
                            0, 1-r(1)-s(1),    0, r(1),    0, s(1) ];
            N2 = [1-r(2)-s(2),           0, r(2),    0, s(2),   0   ;
                            0, 1-r(2)-s(2),    0, r(2),    0, s(2) ];
            N3 = [1-r(3)-s(3),           0, r(3),    0, s(3),   0   ;
                            0, 1-r(3)-s(3),    0, r(3),    0, s(3) ];
                        
            % mass matrix computation (3-point integration)
            Me = w3*(N1.'*N1 + N2.'*N2 + N3.'*N3); % area-normalized mass matrix
            obj.M = obj.area.*repmat(Me,[1,1,nTri]);
            
            % stiffness matrix computation (1-point integration)
            obj.K = obj.area.* (obj.B.' * obj.B);

            %             
%             % strain rate - velocity matrix
%             B = obj.dNdX(:,:,:)
%             
        end
        
    end
    
end

