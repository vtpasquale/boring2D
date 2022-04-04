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
        function [DuUDx,DuUDy] = computeConvectionDerivative(obj,ui,Uj)
            % Compute convection derivatives d(ui*Uj)/dx and d(ui*Uj)/dy at 
            % element center from nodal values of ui and Ui
            %
            % Inputs
            % u = [nNodes,1 double] velocities [u1;u2;...;unNodes]
            % U = [nNodes,1 double] scalar [U1;U2;...;UnNodes]
            
            nTri = size(obj.nodeIDs,1);
            u_e = Matrix3D( permute( ui(obj.nodeIDs.'),[1,3,2]) );
            U_e = Matrix3D( permute( Uj(obj.nodeIDs.'),[1,3,2]) );
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
        function [div_uU1,div_uU2] = computeDivergenceUsingConvectionDerivatives(obj,u,U)
            % Compute divergence of [u]*Ui at element center from nodal
            % values of u and U - using convection derivatives 
            % (this function should be used for testing only)
            %
            % Inputs
            % u = [2*nNodes,1 double] velocities [u1;v1;u2;v2;...;unNodes;vnNodes]
            % U = [2*nNodes,1 double] scalar [U1;V1;U2;V2;...;UnNodes;VnNodes]
            u1 = u(1:2:end,1);
            u2 = u(2:2:end,1);
            U1 = U(1:2:end,1);
            U2 = U(2:2:end,1);
            [du1U1dx,~] = obj.computeConvectionDerivative(u1,U1);
            [~,du2U1dy] = obj.computeConvectionDerivative(u2,U1);
            [du1U2dx,~] = obj.computeConvectionDerivative(u1,U2);
            [~,du2U2dy] = obj.computeConvectionDerivative(u2,U2);
            div_uU1 = du1U1dx + du2U1dy;
            div_uU2 = du1U2dx + du2U2dy;
        end
        function [div_uU1,div_uU2] = computeDivergence(obj,u,U)
            % Compute divergence of [u]*Ui at element center from nodal values of u and U
            %
            % Inputs
            % u = [2*nNodes,1 double] velocities [u1;v1;u2;v2;...;unNodes;vnNodes]
            % U = [2*nNodes,1 double] scalar [U1;V1;U2;V2;...;UnNodes;VnNodes]
            nTri = size(obj.nodeIDs,1);
            
            % element dof
            u_e = Matrix3D( permute( u(obj.gDof.'),[1,3,2]) );
            U_e = Matrix3D( permute( U(obj.gDof.'),[1,3,2]) );
            
            % Gauss integration points & weight factor
            r = [1/3 0 0];
            s = [1/3 0 0];
            [n1,~,~]=obj.shapeFunValsAtIntPoints(r,s);
            N = Matrix3D(repmat(n1,[1,1,nTri]));
            DNDx = obj.dNdx;
            DNDy = obj.dNdy;
            
            % Iteration-constant terms
            % uj [2,6,nTri double]
            % row-component j is summed in the divergence calculation
            uj = N*u_e;
            DujDx = DNDx*u_e;
            DujDy = DNDy*u_e;
            
            % Divergence coefficent term divUi [2,6,nTri double]
            % divergence([u]*Ui) = divuUi*U_e  
            % row-component i corresponds to divergence of independent terms U1, U2
            div_uUi = (DujDx(1,:,:)+DujDy(2,:,:)).*N + uj(1,:,:).*DNDx + uj(2,:,:).*DNDy;
                 
            % compute divergence
            divergenceuUi = div_uUi*U_e;
            div_uU1 = squeeze( divergenceuUi(1,:,:) );
            div_uU2 = squeeze( divergenceuUi(2,:,:) );
        end
        function obj = computeConvectionMatrices(obj,u)
            % Computes convection-related matrices
            nTri = size(obj.nodeIDs,1);
            
            % element dof
            u_e = Matrix3D( permute( u(obj.gDof.'),[1,3,2]) );
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            [n1,n2,n3]=obj.shapeFunValsAtIntPoints(r,s);
            N1 = Matrix3D(repmat(n1,[1,1,nTri]));
            N2 = Matrix3D(repmat(n2,[1,1,nTri]));
            N3 = Matrix3D(repmat(n3,[1,1,nTri]));
            
            % Shape function derivatives are constant inside elements
            DNDx = obj.dNdx;
            DNDy = obj.dNdy;
            
            % Iteration-constant terms
            % uj [2,6,nTri double]
            % row-component j is summed in the divergence calculation
            uj1 = N1*u_e;
            uj2 = N2*u_e;
            uj3 = N3*u_e;
            DujDx = DNDx*u_e;
            DujDy = DNDy*u_e;
            
            % Divergence coefficent terms divUi [2,6,nTri double]
            % divergence([u]*Ui) = divuUi*U_e
            % row-component i corresponds to divergence of independent terms [u]*U1, [u]*U2
            div_uUi1 = (DujDx(1,:,:)+DujDy(2,:,:)).*N1 + uj1(1,:,:).*DNDx + uj1(2,:,:).*DNDy;
            div_uUi2 = (DujDx(1,:,:)+DujDy(2,:,:)).*N2 + uj2(1,:,:).*DNDx + uj2(2,:,:).*DNDy;
            div_uUi3 = (DujDx(1,:,:)+DujDy(2,:,:)).*N3 + uj3(1,:,:).*DNDx + uj3(2,:,:).*DNDy;

            % Convection matrix (3-point integration)
            obj.Cu = w3*obj.area.*(N1.'*div_uUi1 + ...
                                   N2.'*div_uUi2 + ...
                                   N3.'*div_uUi3);            
            
            % Convection stabilization matrix (3-point integration)
            obj.Ku = -0.5*w3*obj.area.*( Matrix3D(div_uUi1).'*div_uUi1 + ...
                                         Matrix3D(div_uUi2).'*div_uUi2 + ...
                                         Matrix3D(div_uUi3).'*div_uUi3 );
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

