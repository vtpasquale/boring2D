clear all; close all; clc

%% Load mesh file
gmf = Gmf(fullfile('naca0012','mesh_NACA0012_inv.mesh'));
gmf.writeVTK('naca0012Divergence.vtk');
tri2D = Tri2D(gmf);
tri2D = computeMassMatrix(tri2D);

%% Manufactured solution
x = gmf.nodes(:,1);
y = gmf.nodes(:,2);
a = manufacturedSolution(x,y);

xc = tri2D.computeCenterFromNodes(x);
yc = tri2D.computeCenterFromNodes(y);
ac = manufacturedSolution(xc,yc);

%% Compute detervatives
% using shape functions
[DuUDx,DuUDy] = tri2D.computeConvectionDerivative(a.u1,a.U1);

% analytic
DuUDxA = ac.du1U1dx;
DuUDxD = ac.du1U1dx - DuUDx;

%% Create global vectors from manufactured components
u = zeros(2*size(gmf.nodes,1),1);
U = zeros(2*size(gmf.nodes,1),1);
u(1:2:end,1) = a.u1;
u(2:2:end,1) = a.u2;
U(1:2:end,1) = a.U1;
U(2:2:end,1) = a.U2;

%% Compute divergence
[divU1b,divU2b] = tri2D.computeDivergence(u,U);

% analytic error
divU1Error = abs(ac.divU1 - divU1b);
divU2Error = abs(ac.divU2 - divU2b);

%% Compute divergence indirectly using derivatives
[divU1c,divU2c] = tri2D.computeDivergenceUsingConvectionDerivatives(u,U);

%% Append VTK file
fid = fopen('naca0012Divergence.vtk','a+');

fprintf(fid,'CELL_DATA %d\n',size(gmf.tri,1));
fprintf(fid,'SCALARS uU float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.u1.*ac.U1);

fprintf(fid,'SCALARS divU1a float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.divU1);
fprintf(fid,'SCALARS divU1b float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU1b);
fprintf(fid,'SCALARS divU1c float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU1c);
fprintf(fid,'SCALARS divU1Error float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.divU1-divU1b);

fprintf(fid,'SCALARS divU2a float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.divU2);
fprintf(fid,'SCALARS divU2b float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU2b);
fprintf(fid,'SCALARS divU2c float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU2c);
fprintf(fid,'SCALARS divU2Error float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.divU2-divU2b);

fclose(fid);
