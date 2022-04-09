classdef Tri2dScalarConvectionDiffusion < Tri2D
    % Class for 2D triangular elements for multi-dimensional scalar convection-diffusion
    
    % Fundamentals of the Finite Element Method for Heat and Mass Transfer, Second Edition
    % P. Nithiarasu, R. W. Lewis, and K. N. Seetharamu
    %
    % Section 7.4.3 - Multi-dimensional scalar convection-diffusion
    
    % Anthony Ricciardi
    % April 2022
    
    properties       
        Mcd % [3,3,nTri double] Convection-diffusion element mass matrix
        Kecd % [3,3,nTri double] Convection-diffusion element diffusion matrix (k normalized)
        Ccd % [3,3,nTri double] Convection-diffusion element convection matrix
        Kscd % [3,3,nTri double] Convection-diffusion element stabilization matrix (timestep normalized)
    end
    methods
        function obj = Tri2dScalarConvectionDiffusion(gmf)
            
            % Use parent class constructor
            obj@Tri2D(gmf);
        end
        
        function obj = computeElementMatrices(obj,u)
            [nu,mu]=size(u);
            if nu ~=2; error('This is ment for constant velocity'); end
            if mu ~=1; error('mu ~=1'); end
            
            nTri = size(obj.nodeIDs,1);
            
            % 3-Point Gauss integration points & weight factor
            r = [2/3 1/6 1/6];
            s = [1/6 1/6 2/3];
            w3 = 1/3;
            
            %% Compute mass matrix
            % % Using Gauss integration
            % Shape function values at integration points
            [n1,n2,n3]=obj.scalarShapeFunValsAtIntPoints(r,s);
            
            % Mass matrix (3-point integration)
            Me = w3*(n1.'*n1 + n2.'*n2 + n3.'*n3); % area-normalized mass matrix
            Mgauss = obj.area.*repmat(Me,[1,1,nTri]);
            
            % % Using area coordinates (NLS Eq. 7.119)
            Marea = (1./12)*obj.area.*repmat([2 1 1; 1 2 1; 1 1 2],[1,1,nTri]);
            
            % Mgauss - Marea = 0 % check
            
            %% Compute convection matrix
            % % Using Gauss integration
            N1 = Matrix3D(repmat(n1,[1,1,nTri]));
            N2 = Matrix3D(repmat(n2,[1,1,nTri]));
            N3 = Matrix3D(repmat(n3,[1,1,nTri]));
            C1 = N1.'*obj.dNpdx + N2.'*obj.dNpdx + N3.'*obj.dNpdx;
            C2 = N1.'*obj.dNpdy + N2.'*obj.dNpdy + N3.'*obj.dNpdy;
            Cgauss = w3.*obj.area.*( u(1)*C1 + u(2)*C2 );
            
            % % Using area coordinates (NLS Eq. 7.120)
            [bm,cm] = obj.computerAreaCoordinateMatrices();
            Carea = (1/6)*( u(1).*bm + u(2).*cm );
            
            % Cgauss - Carea = 0 % check
            
            %% Diffusion matrix (k normalized)
            % % Using Gauss integration (1-point)
            KeGauss = obj.area.*(obj.dNpdx.'*obj.dNpdx + obj.dNpdy.'*obj.dNpdy);
            
            % % Using area coordinates (NLS Eq. 7.122)
            KeArea = (1./(4*obj.area)).*( bm.'.*bm + cm.'.*cm );
            
            % KeGauss - KeArea = 0 % check
            
            %% Stabilization matrix (timestep normalized)
            % % Using Gauss integration (1-point)
            KsGauss = 0.5*obj.area.*(...
                u(1).*(u(1).*obj.dNpdx.'*obj.dNpdx + u(2).*obj.dNpdx.'*obj.dNpdy) + ...
                u(2).*(u(1).*obj.dNpdy.'*obj.dNpdx + u(2).*obj.dNpdy.'*obj.dNpdy) );
            
            % % Using area coordinates (NLS Eq. 7.123)
            KsArea = 0.5*(1./(4*obj.area)).*(...
                u(1).*(u(1).*bm.'.*bm + u(2).*bm.'.*cm ) + ...
                u(2).*(u(1).*cm.'.*bm + u(2).*cm.'.*cm ) );
            
            % KsGauss - KsArea = 0 % check
            
            %% Save to object
            obj.Mcd = Mgauss;
            obj.Kecd = KeGauss;
            obj.Ccd = Cgauss;
            obj.Kscd = KsGauss;
            
        end
        function [M,Ke,C,dtKs] = assembleGlobalMatrices(obj,deltaTimeElement)
            nTri = size(obj.nodeIDs,1);
            
            % timestep scaling
            dt = zeros(1,1,nTri);
            dt(:) = deltaTimeElement;
            dtKscd = dt.*obj.Kscd;
            
            % 3x3 Index management
            [rowIndex3x3,colIndex3x3] = compute3x3Indices(obj);
            
            % Assemble 3x3 element matrices
            M = sparse(rowIndex3x3,colIndex3x3,obj.Mcd(:));
            Ke = sparse(rowIndex3x3,colIndex3x3,obj.Kecd(:));
            C = sparse(rowIndex3x3,colIndex3x3,obj.Ccd(:));
            dtKs = sparse(rowIndex3x3,colIndex3x3,dtKscd(:));
        end
    end
end

