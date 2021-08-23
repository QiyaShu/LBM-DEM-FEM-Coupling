clc;
clear all;
%function T = particleThermalConductionSetup(SizeScaling, Density, HeatCapacity, Conductivity)

[p,e,t] = importMeshGmsh3D('sphere.msh');
% Convert mesh to second order
%[p,e,t,nVnodes,nPnodes,indices] = convertMeshToSecondOrder(p,e,t);

SizeScaling=[5.11328272e-5, 5.11328272e-5, 5.11328272e-5];%radius
Density=3970;
HeatCapacity=765;
Conductivity=[25;25;25];

p(1,:) = p(1,:)*SizeScaling(1); %radius a in x direction
p(2,:) = p(2,:)*SizeScaling(2); %radius b in y direction
p(3,:) = p(3,:)*SizeScaling(3); %radius c in z direction

% Define timestepping parameters
%nSubstep = 1; % number of substeps
%dt = TimeStep/nSubstep;

SurfacePoints = unique(e(1:3,e(4,:)==1),'stable');

% Define material parameters
rho = Density; % density
cp = HeatCapacity; % heat capacity
lambda = Conductivity; % thermal conductivity as a 1-by-nelements vector to be consant across a single element, or as a nnodes-by-1 vector which specifies diffusivity in each node

% Precompute mass matrix and diffusion matrix
M = assembleMassMatrix(p,t,rho*cp);
D = assembleDiffusionMatrix3D(p,t,lambda); 

save('matlabSetup_VaSp.mat','p','e','t','M','D','SurfacePoints');
T = 0.5*ones (size(p,2),1);
save('matlabTResult_0.mat','T');
%end