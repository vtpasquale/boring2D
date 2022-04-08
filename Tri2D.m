classdef Tri2D
    % Class for 2D triangular elements
    
    % Anthony Ricciardi
    % March 2022
    
    properties
        nodeIDs % [nTri,3 uint32] node numbers
        x % [nTri,3 double] x node locations
        y % [nTri,3 double] y node locations
        
        gDof % [nTri,6 uint32] global DOF from element DOF. Assumes two dof per node in node order.

        invJ % [2,2,nTri Matrix3D] inverse of Jacobian matrix (constant inside elements)
        area % [:,:,nTri double] element area
        dNdx % [2,6,nTri double] physical shape function x derivatives for vector assembly (constant inside elements)
        dNdy % [2,6,nTri double] physical shape function y derivatives for vector assembly (constant inside elements)
        dNpdx % [1,3,nTri double] physical shape function x derivatives for scalar assembly (constant inside elements)
        dNpdy % [1,3,nTri double] physical shape function y derivatives for scalar assembly (constant inside elements)
        B    % [3,6,nTri Matrix3D] velocity-strain rate matrix (constant inside elements)
        
        M    % [6,6,nTri double] mass matrix
        Ktau % [6,6,nTri double] deviatoric stress stiffness matrix
        Cu   % [6,6,nTri double] Convection matrix
        Ku   % [6,6,nTri double] Convection stabilization matrix
                
        P % [~,~,nTri double] Pressure-velocity stabalization matrix
        H % [3,3,nTri double] Pressure stabalization matrix
        Mp % [3,3,nTri double] Pressure mass matrix            
        G  % [3,6,nTri double] Pressure-velocity coupling matrix (3-point integration)
        
        minimumHeight % [nTri,1 double] Minimum element height
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
            
            obj.dNpdx = Matrix3D(dNdX(1,:,:));
            obj.dNpdy = Matrix3D(dNdX(2,:,:));
            
            % Deviatoric stress stiffness matrix (1-point integration)
            mu = ones(1,1,nTri); % UPDATE PLACEHOLDER !!!!!!!
            m = [1, 1, 0].';
            I0 = diag([2, 2, 1]);
            stressStrain0 = I0-(2/3)*(m*m.');
            stressStrain = mu.* repmat(stressStrain0,[1,1,nTri]);
            obj.Ktau = obj.area.* (obj.B.'*stressStrain*obj.B);
            
            
            obj = obj.computeElementHeight();
        end
 
        function obj = computeMassMatrix(obj)
            % Compute mass matrix
            nTri = size(obj.nodeIDs,1);
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            [N1,N2,N3]=obj.vectorShapeFunValsAtIntPoints(r,s);
                        
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
            [n1,~,~]=obj.vectorShapeFunValsAtIntPoints(r,s);
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
            % divergence([u]*Ui) = div_uUi*U_e  
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
            [n1,n2,n3]=obj.vectorShapeFunValsAtIntPoints(r,s);
            N1 = Matrix3D(repmat(n1,[1,1,nTri]));
            N2 = Matrix3D(repmat(n2,[1,1,nTri]));
            N3 = Matrix3D(repmat(n3,[1,1,nTri]));
            
            % Shape function derivatives are constant inside elements
            DNDx = obj.dNdx;
            DNDy = obj.dNdy;
            
            % Iteration-constant nonlinear terms
            % uj [2,6,nTri double]
            % row-component j is summed in the divergence calculation
            uj1 = N1*u_e;
            uj2 = N2*u_e;
            uj3 = N3*u_e;
            DujDx = DNDx*u_e;
            DujDy = DNDy*u_e;
            
            % Divergence coefficent terms GradTuN [2,6,nTri double]
            % divergence([u]*Ui) = [GradT*(u*N)] *U_e = GradTuN *U_e
            % row-component i corresponds to divergence of independent terms [u]*U1, [u]*U2
            GradTuN1 = (DujDx(1,:,:)+DujDy(2,:,:)).*N1 + uj1(1,:,:).*DNDx + uj1(2,:,:).*DNDy;
            GradTuN2 = (DujDx(1,:,:)+DujDy(2,:,:)).*N2 + uj2(1,:,:).*DNDx + uj2(2,:,:).*DNDy;
            GradTuN3 = (DujDx(1,:,:)+DujDy(2,:,:)).*N3 + uj3(1,:,:).*DNDx + uj3(2,:,:).*DNDy;

            % Convection matrix (3-point integration)
            obj.Cu = w3*obj.area.*(N1.'*GradTuN1 + ...
                                   N2.'*GradTuN2 + ...
                                   N3.'*GradTuN3);            
            
            % Convection stabilization matrix (3-point integration)
            obj.Ku = -0.5*w3*obj.area.*( Matrix3D(GradTuN1).'*GradTuN1 + ...
                                         Matrix3D(GradTuN2).'*GradTuN2 + ...
                                         Matrix3D(GradTuN3).'*GradTuN3 );
                                     
            % Pressure-velocity stabalization matrix (3-point integration)
            GradNp = zeros(2,3,nTri);
            GradNp(1,:,:) = obj.dNpdx;
            GradNp(2,:,:) = obj.dNpdy;
            obj.P =  w3*obj.area.*(Matrix3D(GradTuN1).'*GradNp + ...
                                   Matrix3D(GradTuN2).'*GradNp + ...
                                   Matrix3D(GradTuN3).'*GradNp );                              
        end
        function obj = computeIncompressibleConvectionMatrices(obj,u)
            % Computes convection-related matrices
            nTri = size(obj.nodeIDs,1);
            
            % element dof
            u_e = Matrix3D( permute( u(obj.gDof.'),[1,3,2]) );
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            [n1,n2,n3]=obj.vectorShapeFunValsAtIntPoints(r,s);
            N1 = Matrix3D(repmat(n1,[1,1,nTri]));
            N2 = Matrix3D(repmat(n2,[1,1,nTri]));
            N3 = Matrix3D(repmat(n3,[1,1,nTri]));
            
            % Shape function derivatives are constant inside elements
            DNDx = obj.dNdx;
            DNDy = obj.dNdy;
            
            % Iteration-constant nonlinear terms
            % uj [2,6,nTri double]
            % row-component j is summed in the divergence calculation
            uj1 = N1*u_e;
            uj2 = N2*u_e;
            uj3 = N3*u_e;
            DujDx = DNDx*u_e;
            DujDy = DNDy*u_e;
            
            % x momentum equation
            udN1 = uj1(1,:,:).*DNDx + uj1(2,:,:).*DNDy;
            udN2 = uj2(1,:,:).*DNDx + uj2(2,:,:).*DNDy;
            udN3 = uj3(1,:,:).*DNDx + uj3(2,:,:).*DNDy;

            % Convection matrix (3-point integration)
            obj.Cu = w3*obj.area.*(N1.'*udN1 + ...
                                   N2.'*udN2 + ...
                                   N3.'*udN3);            
            
            % Convection stabilization matrix (3-point integration)
            obj.Ku =      w3*obj.area.*( Matrix3D(udN1).'*udN1 + ...
                                         Matrix3D(udN1).'*udN2 + ...
                                         Matrix3D(udN1).'*udN3 );
                                     
            % Pressure-velocity stabalization matrix (3-point integration)
            GradNp = zeros(2,3,nTri);
            GradNp(1,:,:) = obj.dNpdx;
            GradNp(2,:,:) = obj.dNpdy;
            obj.P =  w3*obj.area.*(Matrix3D(udN1).'*GradNp + ...
                                   Matrix3D(udN2).'*GradNp + ...
                                   Matrix3D(udN3).'*GradNp );                              
        end
        
        function obj = computeIncompressibleConvectionMatrices2(obj,u)
            
            nTri = size(obj.nodeIDs,1);
            
            b = [obj.y(:,2)-obj.y(:,3),obj.y(:,3)-obj.y(:,1),obj.y(:,1)-obj.y(:,2)];
            c = [obj.x(:,3)-obj.x(:,2),obj.x(:,1)-obj.x(:,3),obj.x(:,2)-obj.x(:,1)];
            u1 = u(1:2:end);
            u2 = u(2:2:end);
            u1e = u1(obj.nodeIDs);
            u2e = u2(obj.nodeIDs);
            usu_ =  sum(u1e,2);
            vsv_ =  sum(u2e,2);
            usu = zeros(1,1,nTri);
            usu(:) = usu_;
            vsv = zeros(1,1,nTri);
            vsv(:) = vsv_;
            
            bm = Matrix3D(zeros(3,3,nTri));
            cm = bm;
            bm(1,1,:) = b(:,1);
            bm(2,1,:) = b(:,1);
            bm(3,1,:) = b(:,1);
            bm(1,2,:) = b(:,2);
            bm(2,2,:) = b(:,2);
            bm(3,2,:) = b(:,2);
            bm(1,3,:) = b(:,3);
            bm(2,3,:) = b(:,3);
            bm(3,3,:) = b(:,3);
            cm(1,1,:) = c(:,1);
            cm(2,1,:) = c(:,1);
            cm(3,1,:) = c(:,1);
            cm(1,2,:) = c(:,2);
            cm(2,2,:) = c(:,2);
            cm(3,2,:) = c(:,2);
            cm(1,3,:) = c(:,3);
            cm(2,3,:) = c(:,3);
            cm(3,3,:) = c(:,3);
            
            u1m = zeros(3,3,nTri);
            u2m = u1m;
            u1m(1,1,:) = u1e(:,1);
            u1m(1,2,:) = u1e(:,1);
            u1m(1,3,:) = u1e(:,1);
            u1m(2,1,:) = u1e(:,2);
            u1m(2,2,:) = u1e(:,2);
            u1m(2,3,:) = u1e(:,2);
            u1m(3,1,:) = u1e(:,3);
            u1m(3,2,:) = u1e(:,3);
            u1m(3,3,:) = u1e(:,3);
            u2m(1,1,:) = u2e(:,1);
            u2m(1,2,:) = u2e(:,1);
            u2m(1,3,:) = u2e(:,1);
            u2m(2,1,:) = u2e(:,2);
            u2m(2,2,:) = u2e(:,2);
            u2m(2,3,:) = u2e(:,2);
            u2m(3,1,:) = u2e(:,3);
            u2m(3,2,:) = u2e(:,3);
            u2m(3,3,:) = u2e(:,3);
            
            % Convection matrix
            % Matches result from computeIncompressibleConvectionMatrices()
            C = (1/24) * ( (usu+u1m).*bm + (vsv+u2m).*cm );
            obj.Cu(1:2:end,1:2:end,:) = C;
            obj.Cu(2:2:end,2:2:end,:) = C; % 

            % Stabilization matrix            
            u1av_  = u1e(:,1).*(u1e(:,1)+usu(:)) + u1e(:,2).*(u1e(:,2)+usu(:)) + u1e(:,3).*(u1e(:,3)+usu(:));
            u12av_ = u1e(:,1).*(u2e(:,1)+vsv(:)) + u1e(:,2).*(u2e(:,2)+vsv(:)) + u1e(:,3).*(u2e(:,3)+vsv(:));
            u2av_  = u2e(:,1).*(u2e(:,1)+vsv(:)) + u2e(:,2).*(u2e(:,2)+vsv(:)) + u2e(:,3).*(u2e(:,3)+vsv(:));
            
            % u21av_ = u2e(:,1).*(u1e(:,1)+usu(:)) + u2e(:,2).*(u1e(:,2)+usu(:)) + u2e(:,3).*(u1e(:,3)+usu(:));
            u1av  = zeros(1,1,nTri); u1av(:)  = u1av_;
            u2av  = zeros(1,1,nTri); u2av(:)  = u2av_;
            u12av = zeros(1,1,nTri); u12av(:) = u12av_;
            
            Ks = (1./(48*obj.area)).*(  u1av.*(bm.'.*bm) + u12av.*((bm.'.*cm)+(cm.'.*bm)) + u2av.*(cm.'.*cm));
            obj.Ku = zeros(6,6,nTri);
            obj.Ku(1:2:end,1:2:end,:) = Ks;
            obj.Ku(2:2:end,2:2:end,:) = Ks;
            
            %  Momentum diffusion
            Km = (1./(4*obj.area)).*(  u1av.*(bm.'.*bm) + u2av.*(cm.'.*cm));
            obj.Ktau = zeros(6,6,nTri);
            obj.Ktau(1:2:end,1:2:end,:) = Km;
            obj.Ktau(2:2:end,2:2:end,:) = Km;
            
            
        end
        
        function obj = computePressureEquationMatrices(obj,betaIn)
            nTri = size(obj.nodeIDs,1);
            
            % Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            % Shape function values at integration points
            [np1,np2,np3]=obj.scalarShapeFunValsAtIntPoints(r,s);
            [nu1,nu2,nu3]=obj.vectorShapeFunValsAtIntPoints(r,s);
            sumNu1 = Matrix3D(repmat(sum(nu1),[1,1,nTri]));
            sumNu2 = Matrix3D(repmat(sum(nu2),[1,1,nTri]));
            sumNu3 = Matrix3D(repmat(sum(nu3),[1,1,nTri]));
            
            % Pressure mass matrix (3-point integration)
            mp = w3*(np1.'*np1 + np2.'*np2 + np3.'*np3); % area-normalized mass matrix
            beta = zeros(1,1,nTri);
            beta(:) = betaIn;
            obj.Mp = (1/beta.^2).*obj.area.*repmat(mp,[1,1,nTri]);
            
            % Shape function derivatives are constant inside elements
            DNpDx = obj.dNpdx;
            DNpDy = obj.dNpdy;
            
            % Pressure stabalization matrix (1-point integration)
            obj.H = obj.area.* (DNpDx.'*DNpDx + DNpDy.'*DNpDy);
            
            % Pressure-velocity coupling matrix (3-point integration)
            obj.G = w3*obj.area.*(DNpDx.'*sumNu1 + DNpDy.'*sumNu1 + ...
                                  DNpDx.'*sumNu2 + DNpDy.'*sumNu2 + ...
                                  DNpDx.'*sumNu3 + DNpDy.'*sumNu3 );
        end
        
        function [Mu,Cu,Ktau,dtKu,Mp,H,G,P] = assembleGlobalMatrices(obj,deltaTimeElement)
            nTri = size(obj.nodeIDs,1);
            
            % 6x6 Index management
            colDof6x6 = zeros(6,6,nTri,'uint32');
            rowDof6x6 = colDof6x6;
            colDof6x6(1,:,:) = permute(obj.gDof,[3,2,1]);
            colDof6x6(2,:,:) = colDof6x6(1,:,:);
            colDof6x6(3,:,:) = colDof6x6(1,:,:);
            colDof6x6(4,:,:) = colDof6x6(1,:,:);
            colDof6x6(5,:,:) = colDof6x6(1,:,:);
            colDof6x6(6,:,:) = colDof6x6(1,:,:);
            colIndex6x6 = double( colDof6x6(:) );
            
            rowDof6x6(:,1,:) = permute(obj.gDof,[2,3,1]);
            rowDof6x6(:,2,:) = rowDof6x6(:,1,:);
            rowDof6x6(:,3,:) = rowDof6x6(:,1,:);
            rowDof6x6(:,4,:) = rowDof6x6(:,1,:);
            rowDof6x6(:,5,:) = rowDof6x6(:,1,:);
            rowDof6x6(:,6,:) = rowDof6x6(:,1,:);
            rowIndex6x6 = double( rowDof6x6(:) );
            
            % Assemble 6x6 element matrices
            Mu = sparse(rowIndex6x6,colIndex6x6,obj.M(:));
            Cu = sparse(rowIndex6x6,colIndex6x6,obj.Cu(:));
            Ktau = sparse(rowIndex6x6,colIndex6x6,obj.Ktau(:));
            
            dt = zeros(1,1,nTri);
            dt(:) = deltaTimeElement;
            dtKu_e = dt.*obj.Ku;
            dtKu = sparse(rowIndex6x6,colIndex6x6,dtKu_e(:));
            
            % 3x3 Index management
            colDof3x3 = zeros(3,3,nTri,'uint32');
            rowDof3x3 = colDof3x3;
            colDof3x3(1,:,:) = permute(obj.nodeIDs,[3,2,1]);
            colDof3x3(2,:,:) = colDof3x3(1,:,:);
            colDof3x3(3,:,:) = colDof3x3(1,:,:);
            colIndex3x3 = double( colDof3x3(:) );
            
            rowDof3x3(:,1,:) = permute(obj.nodeIDs,[2,3,1]);
            rowDof3x3(:,2,:) = rowDof3x3(:,1,:);
            rowDof3x3(:,3,:) = rowDof3x3(:,1,:);
            rowIndex3x3 = double( rowDof3x3(:) );
            
            % Assemble 3x3 element matrices
            Mp = sparse(rowIndex3x3,colIndex3x3,obj.Mp(:));
            H = sparse(rowIndex3x3,colIndex3x3,obj.H(:));
            
            % 3x6 Index management
            colDof3x6 = zeros(3,6,nTri,'uint32');
            rowDof3x6 = colDof3x6;
            colDof3x6(1,:,:) = permute(obj.gDof,[3,2,1]);
            colDof3x6(2,:,:) = colDof3x6(1,:,:);
            colDof3x6(3,:,:) = colDof3x6(1,:,:);
            colIndex3x6 = double( colDof3x6(:) );
            
            rowDof3x6(:,1,:) = permute(obj.nodeIDs,[2,3,1]);
            rowDof3x6(:,2,:) = rowDof3x6(:,1,:);
            rowDof3x6(:,3,:) = rowDof3x6(:,1,:);
            rowDof3x6(:,4,:) = rowDof3x6(:,1,:);
            rowDof3x6(:,5,:) = rowDof3x6(:,1,:);
            rowDof3x6(:,6,:) = rowDof3x6(:,1,:);
            rowIndex3x6 = double( rowDof3x6(:) );
            
            % Assemble 3x6 element matrices
            G = sparse(rowIndex3x6,colIndex3x6,obj.G(:));
            
            % 6x3 Index management
            colDof6x3 = zeros(6,3,nTri,'uint32');
            rowDof6x3 = colDof6x3;
            colDof6x3(1,:,:) = permute(obj.nodeIDs,[3,2,1]);
            colDof6x3(2,:,:) = colDof6x3(1,:,:);
            colDof6x3(3,:,:) = colDof6x3(1,:,:);
            colDof6x3(4,:,:) = colDof6x3(1,:,:);
            colDof6x3(5,:,:) = colDof6x3(1,:,:);
            colDof6x3(6,:,:) = colDof6x3(1,:,:);
            colIndex6x3 = double( colDof6x3(:) );
            
            rowDof6x3(:,1,:) = permute(obj.gDof,[2,3,1]);
            rowDof6x3(:,2,:) = rowDof6x3(:,1,:);
            rowDof6x3(:,3,:) = rowDof6x3(:,1,:);
            rowIndex6x3 = double( rowDof6x3(:) );
            
            % Assemble 6x3 element matrices
            P = sparse(rowIndex6x3,colIndex6x3,obj.P(:));
        end
    end
    methods (Access = private)
        function obj = computeElementHeight(obj)
            % Compute element minimum height for time step calculations
            
            % using geometry
            area_ = obj.area(:);
            length1 = sqrt((obj.x(:,2)-obj.x(:,3)).^2 + (obj.y(:,2)-obj.y(:,3)).^2);
            length2 = sqrt((obj.x(:,3)-obj.x(:,1)).^2 + (obj.y(:,3)-obj.y(:,1)).^2);
            length3 = sqrt((obj.x(:,1)-obj.x(:,2)).^2 + (obj.y(:,1)-obj.y(:,2)).^2);
            height = 2*area_./[length1,length2,length3];
            obj.minimumHeight = min(height,[],2);
            
            % % Equivalent calculation
            % % Using shape function derivatives (cbsflow2d_main.f90 line 233)
            % % CBSFlow2dICEX calculates two equivalent ways - the first result
            % % is simply overwritten by the second result in CBSFlow2dICEX.
            % minimumHeight_ = min(sqrt(obj.dNpdx.^2 + obj.dNpdy.^2).^(-1));
            % minimumHeight2 = minimumHeight_(:);
        end
    end
    methods (Static = true, Access = private)
        function [N1,N2,N3]=vectorShapeFunValsAtIntPoints(r,s)
            % 2D vector shape function values at integration points
            % shapeFunctions = @(r,s) [1-r-s, r, s];
            N1 = [1-r(1)-s(1),           0, r(1),    0, s(1),   0   ;
                            0, 1-r(1)-s(1),    0, r(1),    0, s(1) ];
            N2 = [1-r(2)-s(2),           0, r(2),    0, s(2),   0   ;
                            0, 1-r(2)-s(2),    0, r(2),    0, s(2) ];
            N3 = [1-r(3)-s(3),           0, r(3),    0, s(3),   0   ;
                            0, 1-r(3)-s(3),    0, r(3),    0, s(3) ];
        end
        function [N1,N2,N3]=scalarShapeFunValsAtIntPoints(r,s)
            % Scalar shape function values at integration points
            % shapeFunctions = @(r,s) [1-r-s, r, s];
            N1 = [1-r(1)-s(1), r(1), s(1)];
            N2 = [1-r(2)-s(2), r(2), s(2)];
            N3 = [1-r(3)-s(3), r(3), s(3)];
        end
    end

    
end

