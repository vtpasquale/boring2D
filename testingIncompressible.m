clear all; close all; clc

caseName = fullfile('ldc2d-re400','5000NUcav');


%% Load parameters file and boundary file
par = CbsFlowParamaters([caseName,'.par']);
bco = CbsBoundaryDefinition([caseName,'.bco']);

%% Load mesh file and assemble
gmf = Gmf([caseName,'.plt']);

tri2D = Tri2dIncompressibleLaminarNS(gmf);
edge2D = Edge2D(gmf);
nTri = size(tri2D.nodeIDs,1);
nNodes = size(gmf.nodes,1);

%% Initial conditions
if par.restart==1
    % Load restart file
    var = CbsRestartFile([caseName,'.var'],nNodes);
    u1 = var.u1;
    u2 = var.u2;
    p = var.p;
    T = var.T;
else
    % Initial Conditions
    u1 = par.Ux*ones(size(gmf.nodes,1),1);
    u2 = par.Uy*ones(size(gmf.nodes,1),1);
    p = par.P*ones(size(gmf.nodes,1),1);
    T = par.T*ones(size(gmf.nodes,1),1);
end

%% various logic and checks
% Convection Logic
if par.convectionType==0
    ani = 1/par.Re; % for forced and mixed convection problems
else
    error('Test for natural convection')
    ani = par.Pr;
end

% Supported boundary conditions
if any(~or(bco.flagCode~=500,bco.flagCode~=503))
    error('Only boundary conditions types 500 and 503 are supported');
end

% Local time stepping
if par.beta_opt~=0; error('update for par.beta_opt~=0'); end

%% Boundary Condition Processing
edge = edge2D;
for i = 1:size(bco.flagList,1)
    edge.boundaryID( edge2D.boundaryID==bco.flagList(i) ) = bco.flagCode(i);
end
nodes500_ = edge.nodeIDs(edge.boundaryID==500,:);
nodesWall = unique( nodes500_(:));
nodes503_ = edge.nodeIDs(edge.boundaryID==503,:);
nodesLid = unique( nodes503_(:));

tic
%% Begin time stepping
realTime = 0;
startStep = 1;
for i = startStep:par.nRealTimesteps
    realTime = realTime + par.realTimestepSize;
    
    % Calculate real time
    if par.realTimestepSize < 1e10
        error('Add logic for time-accurate analysis.')
    end
    
    % Begin pseudo time stepping
    pseudoTimeStep = 0;
    minRealTimeStep=par.csafm*par.realTimestepSize;
    re_half=0.5/ani;
    two_inv_re=2*ani;
    
    while pseudoTimeStep < 2000 % par.ntime
        pseudoTimeStep = pseudoTimeStep + 1;
        if mod(pseudoTimeStep,10)==0
            fprintf('Pseudotime Step: %d\n',pseudoTimeStep);
        end
        if mod(pseudoTimeStep,100)==0
            figure(1)
            loglog(1:pseudoTimeStep-2,[du1./du1(1);du2./du2(1);dp./dp(1)].')
            legend('du1','du2','p')
            grid on
            drawnow()
        end
        
        tri2D = tri2D.computeElementMatrices(u1,u2);
        velocity = sqrt(u1.^2 + u2.^2 + 0.1E-15);
        
        % track changes
        if pseudoTimeStep > 1
            du1(pseudoTimeStep-1) = max((u10-u1).^2);
            du2(pseudoTimeStep-1) = max((u20-u2).^2);
            dp(pseudoTimeStep-1) = max((p0-p).^2);
        end
        
        % previous values
        u10=u1;
        u20=u2;
        p0=p;
        T0=T;
        
        % Calculate step size for elements
        % if par.beta_opt~=0; error('update for par.beta_opt~=0'); end
        maxElementVelocity_ = max( velocity( permute(tri2D.nodeIDs,[3,2,1]) ) );
        maxElementVelocity = maxElementVelocity_(:);
        beta = max([par.epsilon*ones(nTri,1),...
                    maxElementVelocity,...
                    two_inv_re./tri2D.minimumHeight],[],2);
        deltaTimeElement = min([par.csafm*tri2D.minimumHeight.^2*re_half,...
                                par.csafm*tri2D.minimumHeight./(maxElementVelocity+beta),...
                                minRealTimeStep*ones(nTri,1)],[],2);
        % Assemble
        if pseudoTimeStep == 1
            [Mdt,Mdtb2,C,K,dtK,dtKs,dtP1,dtP2,G1,G2] = tri2D.assembleGlobalMatrices(deltaTimeElement,beta);
        else
            [Mdt,Mdtb2,C,K,dtK,dtKs,dtP1,dtP2] = tri2D.assembleGlobalMatrices(deltaTimeElement,beta);
        end
        
        % Step 1 - CBSflow consistent
        invDiagMdt = full(sum(Mdt,2)).^-1;
        deltaU1 = -invDiagMdt.*(  C+ani*  K+.5* dtKs)*u1;
        deltaU2 = -invDiagMdt.*(  C+ani*  K+.5* dtKs)*u2;
        u1 = u1 + deltaU1;
        u2 = u2 + deltaU2;
        
        % Apply BC
        u1(nodesLid) = 1.0;
        u2(nodesLid) = 0;
        u1(nodesWall) = 0;
        u2(nodesWall) = 0;
        
        % Step 2 - CBSflow consistent
        invDiagMdtbt2 = full(sum(Mdtb2,2)).^-1; % spot on
        % this matrix seems to be divded by dt when it should not be - or there is a book typo - review this
        rhs1 = - dtK*p; % spot on
        rhs2 = -(G1*u10 + G2*u20); % spot on
        rhs3 = -par.theta1*(G1*(u1-u10) + G2*(u2-u20)); % spot on
        rhs = rhs1+rhs2+rhs3;
        p = p+ invDiagMdtbt2.*rhs;
        
        % Step 3 - CBSflow consistent
        deltaU1ss = invDiagMdt.*( G1.'*p0 - 0.5*dtP1*p0);
        deltaU2ss = invDiagMdt.*( G2.'*p0 - 0.5*dtP2*p0);       
        u1 = u1 + deltaU1ss;
        u2 = u2 + deltaU2ss;
        
        % Apply BC
        u1(nodesLid) = 1.0;
        u2(nodesLid) = 0;
        u1(nodesWall) = 0;
        u2(nodesWall) = 0;
    end
end
timeStepsToc = toc

%% Write to VTK file
gmf.writeVTK([caseName,'.vtk']);
fid = fopen([caseName,'.vtk'],'a+');
fprintf(fid,'POINT_DATA %d\n',nNodes);
fprintf(fid,'SCALARS u1 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',u1);
fprintf(fid,'SCALARS u2 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',u2);
fprintf(fid,'SCALARS p float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',p);
fclose(fid);
