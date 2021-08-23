%function [TSurface,T] = particleThermalComputation(nStep,OutStep,dt,HeatFlux,TlastStep)

load('matlabSetup_VaSp.mat');

nStep = 0;
OutStep = 2;
dt=3.6101e-6;

%T = TlastStep;
fileInput = "matlabTResult_"+nStep+".mat";
load(fileInput);
%fileOutput = "matlabTResult_"+(nStep+1)*dt+"s.mat";
fileOutput = "matlabTResult_"+nStep+".mat";

heatflux = zeros(size(p,2),1);
%heatflux(SurfacePoints) = HeatFlux;
heatflux(SurfacePoints) = 10*ones(size(SurfacePoints,1),1);

% Compute transient part of the right hand side
F = M/dt*T;
    
% Assemble problem matrix
K = M/dt + D;

[K, F] = imposeScalarBoundaryCondition(p,e,K,F,1,'flux',heatflux);% HeatFlux can be row vectors with the number of elements equal to the total number of nodes in the mesh
T = K\F;% M(T_(k+1)-T_k)/dt=-DT_(k+1) therefore (M/dt+D)T_(k+1)=(M/dt)T_k


TSurface=T(SurfacePoints,:);
% Display solution
figure(1);
displayMesh(p,t);
hold on
displaySolution(p,t,T,'Temperature [K]');
figure(2);
displaySurface(p,e,[4,8,13,16],T,'Temperature [K]');%x-y plane
figure(3);
displaySurface(p,e,[3,10,14,20],T,'Temperature [K]');%y-z plane
figure(4);
displaySurface(p,e,[5,7,11,18],T,'Temperature [K]');%x-z plane

if nStep~=0 && mod(nStep,OutStep)==0
    save(fileOutput,'T');
    vtkwrite('particle_Temperature.vtk', 'unstructured_grid',p(1,:),p(2,:),p(3,:), 'scalars','Temperature',T','Precision',8);
end

%end