clear all; close all; clc

caseName = fullfile('ldc2d-re400','5000NUcav');


%% Load parameters file and boundary file
par = CbsFlowParamaters([caseName,'.par']);
bco = CbsBoundaryDefinition([caseName,'.bco']);

%% Load mesh file and assemble
gmf = Gmf([caseName,'.plt']);
gmf.writeVTK([caseName,'.vtk']);

tri2D = Tri2D(gmf);
edge2D = Edge2D(gmf);
nTri = size(tri2D.nodeIDs,1);
nNodes = size(gmf.nodes,1);

%% Initial conditions
u = zeros(2*size(gmf.nodes,1),1);
if par.restart==1
    % Load restart file
    var = CbsRestartFile(fullfile('ldc2d-re400','5000NUcav.var'),nNodes);
    u(1:2:end,1) = var.u1;
    u(2:2:end,1) = var.u2;
    p = var.p;
    T = var.T;
else
    % Initial Conditions
    u(1:2:end,1) = par.Ux;
    u(2:2:end,1) = par.Uy;
    p = par.P*ones(size(gmf.nodes,1),1);
    T = par.T*ones(size(gmf.nodes,1),1);
end
velocity = sqrt(u(1:2:end,1).^2 + u(2:2:end,1).^2 + 0.1E-15);

%% Write intial solution to file
fid = fopen([caseName,'.vtk'],'a+');
fprintf(fid,'POINT_DATA %d\n',nNodes);
fprintf(fid,'SCALARS u1 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',u(1:2:end,1));
fprintf(fid,'SCALARS u2 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',u(2:2:end,1));
fprintf(fid,'SCALARS p float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',p);
fclose(fid);

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

%% Compute mass matrix
tri2D = computeMassMatrix(tri2D);
% 
% 
% tri2D = computePressureEquationMatrices(tri2D,999);
% tri2D = computeConvectionMatrices(tri2D,u);
% [Mu,Cu,Ktau,Ku,Mp,H,G,P] = tri2D.assembleGlobalMatrices();
% Mud = full(sum(Mu,2));% lumped mass matrix

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
    pseudoTimeSteps = 0;
    minRealTimeStep=par.csafm*par.realTimestepSize;
    re_half=0.5/ani;
    two_inv_re=2*ani;
    
%     while pseudoTimeSteps < par.ntime
        % previous values
        u1=u;
        p1=p;
        T1=T;
        
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
        % These two terms are equal for me
        % par.csafm*tri2D.minimumHeight.^2*re_half - (par.csafm*tri2D.minimumHeight./(maxElementVelocity+beta) )
        
        % convective matricies
        tri2D = computeConvectionMatrices(tri2D,u);
%         tri2D = tri2D.computeIncompressibleConvectionMatrices(u);
        tri2D = tri2D.computeIncompressibleConvectionMatrices2(u);
%         C = tri2D.Cu(:,:,1);
%         gDof = tri2D.gDof(1,:);
%         usu = 3*u(1);
%         vsv = 3*u(2);
%         bi = tri2D.y(1,2) - tri2D.y(1,3);
%         ci = tri2D.x(1,3) - tri2D.x(1,2);
%         ((1/24)*(vsv+u(2))*ci) - C(1,1)
        tri2D = computePressureEquationMatrices(tri2D,beta);
        [Mu,Cu,Ktau,dtKu,Mp,H,G,P] = tri2D.assembleGlobalMatrices(deltaTimeElement);
        
%         rhs = -Cu*u; % spot on
        

%           rhs = -.5*dtKu*u; % off, insignificant for test case
%         rhs =  -ani*Ktau*u; % % off, significant for test case
        
        rhs = -Cu*u -ani*Ktau*u -.5*dtKu*u;

end

fprintf(1,'%E\n',rhs(1:8))


