clear all; close all; clc

%% Boundary conditions
rho = 1.22500; % kg/m^3
a = 340.294; % m/s - sound speed
uinf = 0.4*a;
UxBoundary = rho*uinf;
UyBoundary = 0;

% mu = 0.0000181206; % [Pa s] dynamic viscosity
% P0 = 101325; % Pa initial pressure
T = 288.150; % K, initial temperature
R = 287.05; % [J/kg K] air gas constant
P0 = rho*R*T;


dt = 0.001/(uinf + a);
theta1 = 0.75;

%% Load mesh file
gmf = Gmf(fullfile('naca0012','mesh_NACA0012_inv.mesh'));
gmf.writeVTK('naca0012.vtk');


%% Assemble mesh
tri2D = Tri2D(gmf);

%% Initial guess
u = zeros(2*size(gmf.nodes,1),1);
u(1:2:end,1) = uinf;
u(2:2:end,1) = 0;
U = rho*u;

p = P0*ones(size(gmf.nodes,1),1);


%% Assemble Constant Element matricies
tri2D = computeMassMatrix(tri2D);
tri2D = computePressureEquationMatrices(tri2D,a);

%% Global system sets
% U is specified at boundary ID #2
USpecifedNodes = gmf.edges(gmf.edges(:,3)==2,1:2);
USpecifedNodes = unique(USpecifedNodes(:));
UxSpecifedDof = 2*USpecifedNodes - 1;
UySpecifedDof = 2*USpecifedNodes;

% Dof sets
nPdof = size(gmf.nodes,1);
nUdof = 2*nPdof;

% U fixed
s = false(nUdof,1);
s(UxSpecifedDof) = true;
s(UySpecifedDof) = true;
numSdof = sum(s);
Us = spalloc(nUdof,1,numSdof);
Us(UxSpecifedDof) = UxBoundary;
Us(UySpecifedDof) = UyBoundary;

% U free
f = ~s;

% for i = 1:100
    % Update global system matrices
    tri2D = computeConvectionMatrices(tri2D,u);
    [Mu,Cu,Ktau,Ku,Mp,H,G,P] = tri2D.assembleGlobalMatrices();
    
    % Partition matrices
    Muff = Mu(f,f);
    Cuff = Cu(f,f);
    Kuff = Ku(f,f);
    pff  = ( Cu(f,s) - dt*Ku(f,s)) *Us(s);
    Uf = U(f);
    
    % Step 1
    deltaUf = -dt*Muff\(Cuff*Uf - dt*Kuff*Uf + pff);
    deltaU1 = zeros(nUdof,1);
    deltaU1(f) = deltaUf;
    
    % Step 2
    deltaP = dt*Mp\(G*U + .5*G*deltaU1 - dt*theta1*H*p);
    
    % Step 3
    deltaU2 = -dt*Mu\(G.'*p + (dt/2)*P*p);
    
    % Increment values
    U = U + (deltaU1 + deltaU2);
    U(s) = Us(s);
    
    p = p + deltaP;
    rho = (1/(R*T)).*p;
    u(1:2:end) = U(1:2:end)./rho;
    u(2:2:end) = U(2:2:end)./rho;
    
% end

%% Append VTK file
fid = fopen('naca0012.vtk','a+');
fprintf(fid,'POINT_DATA %d\n',size(gmf.nodes,1));

fprintf(fid,'SCALARS rho float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',rho);

fprintf(fid,'SCALARS p float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',p);

fprintf(fid,'SCALARS u1 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',u(1:2:end));
fprintf(fid,'SCALARS u2 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',u(2:2:end));

fprintf(fid,'SCALARS U1 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',U(1:2:end));
fprintf(fid,'SCALARS U2 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',U(2:2:end));

fclose(fid);
