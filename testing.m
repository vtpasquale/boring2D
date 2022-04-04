clear all; close all; clc

%% Load mesh file
gmf = Gmf(fullfile('naca0012','mesh_NACA0012_inv.mesh'));
gmf.writeVTK('naca0012.vtk');

%% Assemble mesh
tri2D = Tri2D(gmf);
tri2D = computeMassMatrix(tri2D);

%% Manufactured solution
x = gmf.nodes(:,1);
y = gmf.nodes(:,2);
a = manufacturedSolution(x,y);

xc = tri2D.computeCenterFromNodes(x);
yc = tri2D.computeCenterFromNodes(y);
ac = manufacturedSolution(xc,yc);

% Create global vectors from manufactured components
u = zeros(2*size(gmf.nodes,1),1);
U = zeros(2*size(gmf.nodes,1),1);
u(1:2:end,1) = a.u1;
u(2:2:end,1) = a.u2;
U(1:2:end,1) = a.U1;
U(2:2:end,1) = a.U2;

%% Assemble matricies
tri2D = computeConvectionMatrices(tri2D,u);


% %% Append VTK file
% fid = fopen('naca0012.vtk','a+');
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

% fprintf(fid,'CELL_DATA %d\n',size(gmf.tri,1));
% fprintf(fid,'SCALARS uU float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',ac.u1.*ac.U1);
% 
% fprintf(fid,'SCALARS divU1a float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',ac.divU1);
% fprintf(fid,'SCALARS divU1b float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',divU1b);
% fprintf(fid,'SCALARS divU1c float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',divU1c);
% fprintf(fid,'SCALARS divU1Error float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',ac.divU1-divU1b);
% 
% fprintf(fid,'SCALARS divU2a float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',ac.divU2);
% fprintf(fid,'SCALARS divU2b float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',divU2b);
% fprintf(fid,'SCALARS divU2c float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',divU2c);
% fprintf(fid,'SCALARS divU2Error float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',ac.divU2-divU2b);

% fclose(fid);
