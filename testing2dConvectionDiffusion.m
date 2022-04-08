clear all; close all; clc

caseName = fullfile('ldc2d-re400','5000NUcav');


%% Load mesh file and assemble
gmf = Gmf([caseName,'.plt']);

tri2D = Tri2D(gmf);
edge2D = Edge2D(gmf);
nTri = size(tri2D.nodeIDs,1);
nNodes = size(gmf.nodes,1);

%% Assemble
u = [1;1];
tri2D = tri2D.compute2dConvectionDiffusionMatrices(u);
deltaTimeElement = tri2D.minimumHeight./(sqrt(sum(u.^2)));
[M,Ke,C,dtKs] = tri2D.assemble2dConvectionDiffusionMatrices(deltaTimeElement);

%% Boundary Conditions
leftSide = gmf.nodes(:,1)==0;
bottom = gmf.nodes(:,2)==0;
leftBottom = and(leftSide,gmf.nodes(:,2)<=0.2);
leftTop = and(leftSide,~leftBottom);
phi0 = or(bottom,leftBottom);
phi1 = leftTop;
if any(and(phi0,phi1)); error('Not exclusive'); end
s = or(phi0,phi1);
f = ~s;

phis = zeros(nNodes,1);
phis(phi0) = 0;
phis(phi1) = 1;

%% Static solution
A = [C + dtKs];
F = zeros(nNodes,1);
Aff = A(f,f);
Ff = F(f) - A(f,s)*phis(s);
phif = Aff\Ff;
phi = zeros(nNodes,1);
phi(f) = phif;
phi(s) = phis(s);


%% Write to VTK file
gmf.writeVTK([caseName,'.vtk']);
fid = fopen([caseName,'.vtk'],'a+');
fprintf(fid,'POINT_DATA %d\n',nNodes);
fprintf(fid,'SCALARS phi float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',phi);
fclose(fid);