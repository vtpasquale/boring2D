classdef Tri2dIncompressibleLaminarNS < Tri2D
    % Class for 2D triangular elements for incompressible laminar Navier Stokes equations.
    
    % Fundamentals of the Finite Element Method for Heat and Mass Transfer, Second Edition
    % P. Nithiarasu, R. W. Lewis, and K. N. Seetharamu
    %
    % Section 7.6
    
    % Anthony Ricciardi
    % April 2022
    
    properties
        C % [3,3,nTri double] Convection matrix
        Ks % [3,3,nTri double] Stabilization matrix (timestep normalized)
        M % [3,3,nTri double] Mass matrix
        K % [3,3,nTri double] Diffusion stiffness matrix (inverse-density-normalized pressure diffusion, nu normalized momentum diffusion, k normalized heat diffusion)
        G1 % [3,3,nTri double] Pressure-u1 gradient matrix
        G2 % [3,3,nTri double] Pressure-u2 gradient matrix
        P1 % [3,3,nTri double] Pressure-u1 stabilization matrix
        P2 % [3,3,nTri double] Pressure-u2 stabilization matrix
    end
    methods
        function obj = Tri2dIncompressibleLaminarNS(gmf)
            
            % Use parent class constructor
            obj@Tri2D(gmf);
        end
        
        function obj = computeElementMatrices(obj,u1,u2)
            [nu1,mu1]=size(u1);
            [nu2,mu2]=size(u2);
            if mu1~=1; error('mu1~=1'); end
            if mu2~=1; error('mu2~=1'); end
            if nu1~=nu2; error('nu1~=nu2'); end
            % obj = computeElementMatricesAreaCoordinates(obj,u1,u2);
            obj = computeElementMatricesGaussQuad(obj,u1,u2);
            %             % Compare results
            %             tic
%                         gaussInt = computeElementMatricesGaussQuad(obj,u1,u2);
            %             tocGauss = toc
            %             tic
            %             areaInt = computeElementMatricesAreaCoordinates(obj,u1,u2);
            %             tocArea = toc
            %             max(abs(gaussInt.C(:) - areaInt.C(:)))
            %             max(abs(gaussInt.Ks(:) - areaInt.Ks(:)))
            %             max(abs(gaussInt.M(:) - areaInt.M(:)))
            %             max(abs(gaussInt.K(:) - areaInt.K(:)))
            %             max(abs(gaussInt.G1(:) - areaInt.G1(:)))
            %             max(abs(gaussInt.G2(:) - areaInt.G2(:)))
            
        end
        function [Mdt,dtC,dtK,dt2Ks,dtG1,dtG2,...
                  M,    C,  K, dtKs,  G1,  G2, ...
                  Mdtb2,G1dt,G2dt,dtP1,dtP2] = assembleGlobalMatrices(obj,deltaTimeElement,beta_)
            nTri = size(obj.nodeIDs,1);
            
            % timestep scaling
            dt = zeros(1,1,nTri);
            dt(:) = deltaTimeElement;
            
            beta = zeros(1,1,nTri);
            beta(:) = beta_;
            
            Mdt_   = (1./dt).*obj.M;
            dtC_   = dt.*obj.C;
            dtK_   = dt.*obj.K;
            dtKs_ = dt.*obj.Ks;
            dt2Ks_ = (dt.^2).*obj.Ks;
            dtG1_  = dt.*obj.G1;
            dtG2_  = dt.*obj.G1;
            dtP1_  = dt.*obj.P1;
            dtP2_  = dt.*obj.P2;
            
            Mdtb2_ =(1./(dt.*beta.^2)).*obj.M;
            G1dt_  = (1./dt).*obj.G1;
            G2dt_  = (1./dt).*obj.G1;
            
            % 3x3 Index management
            [rowIndex3x3,colIndex3x3] = compute3x3Indices(obj);
            
            % Assemble 3x3 element matrices
            Mdt = sparse(rowIndex3x3,colIndex3x3,Mdt_(:));
            dtC = sparse(rowIndex3x3,colIndex3x3,dtC_(:));
            dtK = sparse(rowIndex3x3,colIndex3x3,dtK_(:));
            dt2Ks = sparse(rowIndex3x3,colIndex3x3,dt2Ks_(:));
            dtG1 = sparse(rowIndex3x3,colIndex3x3,dtG1_(:));
            dtG2 = sparse(rowIndex3x3,colIndex3x3,dtG2_(:));
            dtP1 = sparse(rowIndex3x3,colIndex3x3,dtP1_(:));
            dtP2 = sparse(rowIndex3x3,colIndex3x3,dtP2_(:));
            
            M = sparse(rowIndex3x3,colIndex3x3,obj.M(:));
            C = sparse(rowIndex3x3,colIndex3x3,obj.C(:));
            K = sparse(rowIndex3x3,colIndex3x3,obj.K(:));
            dtKs = sparse(rowIndex3x3,colIndex3x3,dtKs_(:));
            G1 = sparse(rowIndex3x3,colIndex3x3,obj.G1(:));
            G2 = sparse(rowIndex3x3,colIndex3x3,obj.G2(:));
            
            Mdtb2 = sparse(rowIndex3x3,colIndex3x3,Mdtb2_(:));
            G1dt = sparse(rowIndex3x3,colIndex3x3,G1dt_(:));
            G2dt = sparse(rowIndex3x3,colIndex3x3,G2dt_(:));

        end
    end
    methods (Access = private)
        function obj = computeElementMatricesGaussQuad(obj,u1,u2)
            % Compute element matrices using Gauss quadrature integration.
            % Gauss quadrature method is easier to relate to the original
            % equations and easier to generalize to more element types. It
            % is probably slower than area coordinate integration.
            
            nTri = size(obj.nodeIDs,1);
            [~,mu1]=size(u1);
            [~,mu2]=size(u2);
            if mu1 ~=1; error('mu1 ~=1'); end
            if mu2 ~=1; error('mu2 ~=1'); end
            
            % Element velocities
            u1e = u1(obj.nodeIDs);
            u2e = u2(obj.nodeIDs);
            
            % % For Gauss integration
            [r,s,w3]=obj.gaussPointsAndWeights();
            [n1,n2,n3]=obj.scalarShapeFunValsAtIntPoints(r,s);
            N1 = Matrix3D(repmat(n1,[1,1,nTri]));
            N2 = Matrix3D(repmat(n2,[1,1,nTri]));
            N3 = Matrix3D(repmat(n3,[1,1,nTri]));
            
            % Gauss point velocities
            u1em = zeros(3,1,nTri);
            u1em(:) = u1e.';
            u1g1 = N1*u1em;
            u1g2 = N2*u1em;
            u1g3 = N3*u1em;
            u2em = zeros(3,1,nTri);
            u2em(:) = u2e.';
            u2g1 = N1*u2em;
            u2g2 = N2*u2em;
            u2g3 = N3*u2em;
            
            % Convection matrix  (3-point integration)
            obj.C  = w3*obj.area.* ( N1.'*(u1g1.*obj.dNpdx + u2g1.*obj.dNpdy) + ...
                N2.'*(u1g2.*obj.dNpdx + u2g2.*obj.dNpdy) + ...
                N3.'*(u1g3.*obj.dNpdx + u2g3.*obj.dNpdy) );
            
            
            % Stabilization matrix (3-point)
            obj.Ks = w3*obj.area.*(...
                u1g1.*(u1g1.*obj.dNpdx.'*obj.dNpdx + u2g1.*obj.dNpdx.'*obj.dNpdy) + ...
                u2g1.*(u1g1.*obj.dNpdy.'*obj.dNpdx + u2g1.*obj.dNpdy.'*obj.dNpdy) + ...
                u1g2.*(u1g2.*obj.dNpdx.'*obj.dNpdx + u2g2.*obj.dNpdx.'*obj.dNpdy) + ...
                u2g2.*(u1g2.*obj.dNpdy.'*obj.dNpdx + u2g2.*obj.dNpdy.'*obj.dNpdy) + ...
                u1g3.*(u1g3.*obj.dNpdx.'*obj.dNpdx + u2g3.*obj.dNpdx.'*obj.dNpdy) + ...
                u2g3.*(u1g3.*obj.dNpdy.'*obj.dNpdx + u2g3.*obj.dNpdy.'*obj.dNpdy) );
            
            % Mass matrix (3-point integration)
            Me = w3*(n1.'*n1 + n2.'*n2 + n3.'*n3); % area-normalized
            obj.M = obj.area.*repmat(Me,[1,1,nTri]);
            
            % Diffusion matrix (1-point)
            % + nu-normalized momentum diffusion
            % + k-normalized heat diffusion
            % + inverse-density-normalized pressure diffusion
            obj.K = obj.area.*(obj.dNpdx.'*obj.dNpdx + obj.dNpdy.'*obj.dNpdy);
            
            %% Pressure gradient matrices
            % % Using Gauss integration (1-point) [Ni = Nj = Nk = w3]
            obj.G1  = w3.*obj.area.*[obj.dNpdx;
                obj.dNpdx;
                obj.dNpdx];
            obj.G2  = w3.*obj.area.*[obj.dNpdy;
                obj.dNpdy;
                obj.dNpdy];
            
            %% Pressure Stabilization matrix (3-point integration)
            obj.P1  = w3*obj.area.*(...
                      (u1g1.*obj.dNpdx.'*obj.dNpdx + u2g1.*obj.dNpdy.'*obj.dNpdx) + ...
                      (u1g2.*obj.dNpdx.'*obj.dNpdx + u2g2.*obj.dNpdy.'*obj.dNpdx) + ...
                      (u1g3.*obj.dNpdx.'*obj.dNpdx + u2g3.*obj.dNpdy.'*obj.dNpdx) );
            obj.P2  = w3*obj.area.*(...
                      (u1g1.*obj.dNpdx.'*obj.dNpdy + u2g1.*obj.dNpdy.'*obj.dNpdy) + ...
                      (u1g2.*obj.dNpdx.'*obj.dNpdy + u2g2.*obj.dNpdy.'*obj.dNpdy) + ...
                      (u1g3.*obj.dNpdx.'*obj.dNpdy + u2g3.*obj.dNpdy.'*obj.dNpdy) );
            
        end
        function obj = computeElementMatricesAreaCoordinates(obj,u1,u2)
            % Compute element matrices using Area coordinates
            nTri = size(obj.nodeIDs,1);
            
            % Element velocities
            u1e = u1(obj.nodeIDs);
            u2e = u2(obj.nodeIDs);
            
            % % For area coordinates integration
            [bm,cm] = obj.computerAreaCoordinateMatrices();
            [u1m,u2m] = obj.computeElementVelocityMatrices(u1e,u2e);
            usu_ =  sum(u1e,2);
            vsv_ =  sum(u2e,2);
            usu = zeros(1,1,nTri);
            usu(:) = usu_;
            vsv = zeros(1,1,nTri);
            vsv(:) = vsv_;
            
            % Convection matrix (NLS Eq. 7.153)
            obj.C = (1/24) * ( (usu+u1m).*bm + (vsv+u2m).*cm );
            
            % Stabilization matrix (NLS Eq. 7.157)
            u1av_  = u1e(:,1).*(u1e(:,1)+usu(:)) + u1e(:,2).*(u1e(:,2)+usu(:)) + u1e(:,3).*(u1e(:,3)+usu(:));
            u12av_ = u1e(:,1).*(u2e(:,1)+vsv(:)) + u1e(:,2).*(u2e(:,2)+vsv(:)) + u1e(:,3).*(u2e(:,3)+vsv(:));
            u2av_  = u2e(:,1).*(u2e(:,1)+vsv(:)) + u2e(:,2).*(u2e(:,2)+vsv(:)) + u2e(:,3).*(u2e(:,3)+vsv(:));
            % u21av_ = u2e(:,1).*(u1e(:,1)+usu(:)) + u2e(:,2).*(u1e(:,2)+usu(:)) + u2e(:,3).*(u1e(:,3)+usu(:));
            u1av  = zeros(1,1,nTri); u1av(:)  = u1av_;
            u2av  = zeros(1,1,nTri); u2av(:)  = u2av_;
            u12av = zeros(1,1,nTri); u12av(:) = u12av_;
            obj.Ks = (1./(48*obj.area)).*(  u1av.*(bm.'.*bm) + u12av.*((bm.'.*cm)+(cm.'.*bm)) + u2av.*(cm.'.*cm));
            
            % Compute mass matrix (NLS Eq. 7.152)
            obj.M = (1./12)*obj.area.*repmat([2 1 1; 1 2 1; 1 1 2],[1,1,nTri]);
            
            % Diffusion matrix
            % + nu-normalized momentum diffusion (NLS Eq. 7.155)
            % + k-normalized heat diffusion (NLS Eq. 7.156)
            % + inverse-density-normalized pressure diffusion (NLS Eq. 7.158)
            obj.K = (1./(4*obj.area)).*( bm.'.*bm + cm.'.*cm );
            
            % Pressure gradient matrices (NLS Eq. 7.159-160)
            obj.G1 = (1/6).* bm;
            obj.G2 = (1/6).* cm;
            
        end
    end
end

