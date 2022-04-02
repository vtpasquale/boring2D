clear all; close all; clc

%% Load mesh file
gmf = Gmf(fullfile('naca0012','mesh_NACA0012_inv.mesh'));
gmf.writeVTK('naca0012.vtk');



tri2D = Tri2D(gmf);
tri2D = computeMassMatrix(tri2D);

u = rand(2*size(tri2D.nodeIDs,1),1);
% tri2D = computeConvectionMatrices(tri2D,u);


%% Assumed solution and analytic derivatives
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

%% Divergence from derivatives
[du1U1dx,~] = tri2D.computeConvectionDerivative(a.u1,a.U1);
[~,du2U1dy] = tri2D.computeConvectionDerivative(a.u2,a.U1);

divU1alt = du1U1dx + du2U1dy;

%% Compute divergence
% using shape functions
[divU1,divU2] = tri2D.computeDivergence([a.u1,a.u2],[a.U1,a.U2]);

% analytic error
divU1Error = abs(ac.divU1 - divU1);
divU2Error = abs(ac.divU2 - divU2);

%% Append VTK file
fid = fopen('naca0012.vtk','a+');
% fprintf(fid,'POINT_DATA %d\n',size(gmf.nodes,1));
% fprintf(fid,'SCALARS uUn float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',u1.*U1);
% fprintf(fid,'SCALARS du1U1dx float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',du1U1dx);
% fprintf(fid,'SCALARS du1U1dy float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',du1U1dy);

fprintf(fid,'CELL_DATA %d\n',size(gmf.tri,1));
fprintf(fid,'SCALARS uUc float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.u1.*ac.U1);

fprintf(fid,'SCALARS divU1 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU1);
fprintf(fid,'SCALARS divU1A float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.divU1);
fprintf(fid,'SCALARS divU1Error float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU1Error);
fprintf(fid,'SCALARS divU1alt float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU1alt);

fprintf(fid,'SCALARS divU2 float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU2);
fprintf(fid,'SCALARS divU2A float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',ac.divU2);
fprintf(fid,'SCALARS divU2Error float\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',divU2Error);

% fprintf(fid,'SCALARS du1U1dy float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',DuUDy);

fclose(fid);
