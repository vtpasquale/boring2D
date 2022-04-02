classdef Tri2D
    % Class for 2D triangular elements
    
    % Anthony Ricciardi
    % March 2022
    
    properties
        nodeIDs % [3,nTri uint32] node numbers
        x % [3,nTri double] x node locations
        y % [3,nTri double] y node locations
        
        gDof % [6,nTri uint32] global DOF from element DOF. Assumes two dof per node in node order.

        invJ % [2,2,nTri Matrix3D] inverse of Jacobian matrix (constant inside elements)
        area % [:,:,nTri double] element area
        dNdx % [2,6,nTri double] physical shape function x derivatives for element assembly (constant inside elements)
        dNdy % [2,6,nTri double] physical shape function y derivatives for element assembly (constant inside elements)
        B    % [3,6,nTri Matrix3D] velocity-strain rate matrix (constant inside elements)
        
        M    % [6,6,nTri double] mass matrix
        Ktau % [6,6,nTri double] deviatoric stress stiffness matrix
        Cu   % [6,6,nTri double] Convection matrix
        Ku   % [6,6,nTri double] Convection stabilization matrix
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
            
            % global DOF index (assumes two dof per node in node order)
            obj.gDof =zeros(nTri,6,'uint32');
            obj.gDof(:,[1,3,5]) = 2*obj.nodeIDs - 1;
            obj.gDof(:,[2,4,6]) = 2*obj.nodeIDs;
            
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
            
            % Physical shape function derivative matrix
            dNdX = invJ*dNXiNd;
            
%             % Alternate B computation
%             B2 = zeros(3,6,nTri);
%             B2(1,1,:) = dNdX(1,1,:);
%             B2(1,3,:) = dNdX(1,2,:);
%             B2(1,5,:) = dNdX(1,3,:);
%             B2(2,2,:) = dNdX(2,1,:);
%             B2(2,4,:) = dNdX(2,2,:);
%             B2(2,6,:) = dNdX(2,3,:);
%             B2(3,1,:) = dNdX(2,1,:);
%             B2(3,2,:) = dNdX(1,1,:);
%             B2(3,3,:) = dNdX(2,2,:);
%             B2(3,4,:) = dNdX(1,2,:);
%             B2(3,5,:) = dNdX(2,3,:);
%             B2(3,6,:) = dNdX(1,3,:);
%             B-B2 % check == small
            
            % Physical shape function detervative matrix for element assembly
            obj.dNdx = Matrix3D(zeros(2,6,nTri));
            obj.dNdy = obj.dNdx;
            obj.dNdx(1,[1,3,5],:)=dNdX(1,:,:);
            obj.dNdx(2,[2,4,6],:)=dNdX(1,:,:);
            obj.dNdy(1,[1,3,5],:)=dNdX(2,:,:);
            obj.dNdy(2,[2,4,6],:)=dNdX(2,:,:);
            
            % Deviatoric stress stiffness matrix (1-point integration)
            mu = ones(1,1,nTri); % UPDATE PLACEHOLDER !!!!!!!
            m = [1, 1, 0].';
            I0 = diag([2, 2, 1]);
            stressStrain0 = I0-(2/3)*(m*m.');
            stressStrain = mu.* repmat(stressStrain0,[1,1,nTri]);
            obj.Ktau = obj.area.* (obj.B.'*stressStrain*obj.B);
        end
        
        function obj = computeMassMatrix(obj)
            % Compute mass matrix
            nTri = size(obj.nodeIDs,1);
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            [N1,N2,N3]=obj.shapeFunValsAtIntPoints(r,s);
                        
            % Mass matrix (3-point integration)
            Me = w3*(N1.'*N1 + N2.'*N2 + N3.'*N3); % area-normalized mass matrix
            obj.M = obj.area.*repmat(Me,[1,1,nTri]);
        end
        function u_c = computeCenterFromNodes(obj,u)
            % interpolated value at element center from global scalar
            nTri = size(obj.nodeIDs,1);
            u_e = Matrix3D( permute( u(obj.nodeIDs.'),[1,3,2]) );
            r = 1/3; s = 1/3;
            n = [1-r-s, r, s];
            N = Matrix3D(repmat(n,[1,1,nTri]));
            u_c0 = N*u_e;
            u_c = u_c0(:);
        end
        function [DuUDx,DuUDy] = computeConvectionDerivative(obj,u,U)
            % compute convection derivative at element center from global
            % scalar values of u and U
            nTri = size(obj.nodeIDs,1);
            u_e = Matrix3D( permute( u(obj.nodeIDs.'),[1,3,2]) );
            U_e = Matrix3D( permute( U(obj.nodeIDs.'),[1,3,2]) );
            r = 1/3; s = 1/3;
            n = [1-r-s, r, s];
            N = Matrix3D(repmat(n,[1,1,nTri]));
            DNDx = obj.dNdx(1,[1,3,5],:);
            DNDy = obj.dNdy(1,[1,3,5],:);
            
            Nu = N*u_e;
            DNuDx = DNDx*u_e;
            DNuDy = DNDy*u_e;
            
            % chain rule
            DNuNDx = DNuDx*N + Nu*DNDx;
            DNuNDy = DNuDy*N + Nu*DNDy;
            
            DuUDx0 = DNuNDx*U_e;
            DuUDy0 = DNuNDy*U_e;
            DuUDx  = DuUDx0(:);
            DuUDy  = DuUDy0(:);
        end
        function [divU1,divU2] = computeDivergence(obj,uIn,UIn)
            % compute convection derivative at element center from global
            % scalar values of u and U
            nTri = size(obj.nodeIDs,1);
            
            % Convert to global coordinates
            nG = 2*size(uIn,1);
            u = zeros(nG,1);
            u(1:2:end) = uIn(:,1);
            u(2:2:end) = uIn(:,2);
            U = zeros(nG,1);
            U(1:2:end) = UIn(:,1);
            U(2:2:end) = UIn(:,2);
            
            % element dof
            uTilde = Matrix3D( permute( u(obj.gDof.'),[1,3,2]) );
            UTilde = Matrix3D( permute( U(obj.gDof.'),[1,3,2]) );
            
            % Gauss integration points & weight factor
            r = [1/3 0 0];
            s = [1/3 0 0];
            [n1,n2,n3]=obj.shapeFunValsAtIntPoints(r,s);
            N = Matrix3D(repmat(n1,[1,1,nTri]));
            
            DNDx = obj.dNdx;
            DNDy = obj.dNdy;
            
            uj = N*uTilde;
            DujDx = DNDx*uTilde;
            DujDy = DNDy*uTilde;
            
            % chain rule
            Du1UiDx = DujDx(1,:,:).*N + uj(1,:,:).*DNDx;
            Du2UiDy = DujDy(2,:,:).*N + uj(2,:,:).*DNDy;
            
            % compute divergence
            divU = (Du1UiDx + Du2UiDy)*UTilde;
            divU1 = squeeze( divU(1,:,:) );
            divU2 = squeeze( divU(2,:,:) );
        end
        function obj = computeConvectionMatrices(obj,u)
            % Computes convection-related matrices
            nTri = size(obj.nodeIDs,1);
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            [n1,n2,n3]=obj.shapeFunValsAtIntPoints(r,s);
            N1 = Matrix3D(repmat(n1,[1,1,nTri]));
            N2 = Matrix3D(repmat(n2,[1,1,nTri]));
            N3 = Matrix3D(repmat(n3,[1,1,nTri]));
            
            keyboard
            % update this section

            % Convection matrix (3-point integration)
            uTilde = Matrix3D( permute( u(obj.gDof.'),[1,3,2]) );
            DNu = obj.dNdX(1,:,:)*uTilde + obj.dNdX(2,:,:)*uTilde;
            DuN1 = DNu.*N1 + N1([1,1],:,:)*uTilde*obj.dNdX(1,:,:) + N1([2,2],:,:)*uTilde*obj.dNdX(2,:,:);
            DuN2 = DNu.*N2 + N2([1,1],:,:)*uTilde*obj.dNdX(1,:,:) + N2([2,2],:,:)*uTilde*obj.dNdX(2,:,:);
            DuN3 = DNu.*N3 + N3([1,1],:,:)*uTilde*obj.dNdX(1,:,:) + N3([2,2],:,:)*uTilde*obj.dNdX(2,:,:);
            obj.Cu = w3*obj.area.*(N1.'*DuN1 + N2.'*DuN2 + N3.'*DuN3);            
            
            % Convection stabilization matrix (3-point integration)
            obj.Ku = -0.5*w3*obj.area.*( Matrix3D(DuN1).'*DuN1 + ...
                                         Matrix3D(DuN2).'*DuN2 + ...
                                         Matrix3D(DuN3).'*DuN3 );
        end
        
    end
    methods (Static = true, Access = private)
        function [N1,N2,N3]=shapeFunValsAtIntPoints(r,s)
            % Shape function values at integration points
            % shapeFunctions = @(r,s) [1-r-s, r, s];
            N1 = [1-r(1)-s(1),           0, r(1),    0, s(1),   0   ;
                            0, 1-r(1)-s(1),    0, r(1),    0, s(1) ];
            N2 = [1-r(2)-s(2),           0, r(2),    0, s(2),   0   ;
                            0, 1-r(2)-s(2),    0, r(2),    0, s(2) ];
            N3 = [1-r(3)-s(3),           0, r(3),    0, s(3),   0   ;
                            0, 1-r(3)-s(3),    0, r(3),    0, s(3) ];
        end
    end

    
end

