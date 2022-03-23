clear all; close all; clc

%% Load mesh file
gmf = Gmf(fullfile('naca0012','mesh_NACA0012_inv.mesh'));
gmf.writeVTK('naca0012.vtk');

tic
tri2D = Tri2D(gmf);
tri2D = computeMass(tri2D);
toc