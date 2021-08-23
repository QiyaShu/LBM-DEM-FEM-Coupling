function [D, F] = assembleDiffusionMatrix3D(p, t, nu)

% assembleDiffusionMatrix3D - Assemble global diffusion matrix in 3-D.
%
% This QuickerSim CFD Toolbox function assembles global diffusion matrix.
%
% [D, F] = assembleDiffusionMatrix3D(p, t, nu);
%
% Input arguments:
% p - array of nodal coordinates (see help of the importMeshGmsh function
%     for details).
% t - array of finite elements (see help of the importMeshGmsh function for
%     details).
% nu - diffusivity value. This parameter can be entered in one of the
%      following ways:
%      1. As a constant scalar value and then is used in all mesh elements.
%      2. As a 1-by-nelements vector where nelements is the total
%         number of elements in the mesh. In this case nu is assumed
%         to be consant across a single element only and enables
%         implementation of a nonlinear diffusion problem.
%      3. As a nnodes-by-1 vector which specifies diffusivity in each node.
%
% Output arguments:
% D - assembled, global diffusion matrix with no boundary conditions
%     included.
% F - right-hand side vector of proper size. Because this function does not
%     account for any boundary conditions or source terms, this is always
%     an all zero vector.
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLENAVIERSTOKESMATRIX2D, ASSEMBLESCALARCONVECTIONMATRIX2D,
% ASSEMBLESCALARSOURCETERM2D, ASSEMBLESTOKESMATRIX2D, IMPORTMESHGMSH.

% make sure nu is a row vector
if size(nu, 1) > 3 || (size(nu, 2) > 1 && numel(nu)<=3)
    nu = nu';
end

nelements = size(t,2);
nnodes = size(p,2);
const_nu_flag = false;
element_nu_flag = false;

if length(nu(1,:)) == nnodes
    nu = moveDataToElements(t, nu);
end

% Choose quadrature based on mesh and viscosity interpolation orders
if size(t, 1) == 5
    
    nodes_per_element = 4;
    viscosity_nodes_per_element = 4;
    
    % Quadrature points
    qPoints = [1/4, 1/4, 1/4];
    weights = 1/6;
    
    nqpoints = size(qPoints, 1);
    Nksi = zeros(nodes_per_element, nqpoints);
    Neta = zeros(nodes_per_element, nqpoints);
    Nyip = zeros(nodes_per_element, nqpoints);
    
    % Compute derivatives of all linear shape functions in all qPoints
    Nksi(1, :) = -1*ones(1, nqpoints);
    Nksi(2, :) = 1*ones(1, nqpoints);
    Nksi(3, :) = 0*ones(1, nqpoints);
    Nksi(4, :) = 0*ones(1, nqpoints);
    
    Neta(1, :) = -1*ones(1, nqpoints);
    Neta(2, :) = 0*ones(1, nqpoints);
    Neta(3, :) = 1*ones(1, nqpoints);
    Neta(4, :) = 0*ones(1, nqpoints);
    
    Nyip(1, :) = -1*ones(1, nqpoints);
    Nyip(2, :) = 0*ones(1, nqpoints);
    Nyip(3, :) = 0*ones(1, nqpoints);
    Nyip(4, :) = 1*ones(1, nqpoints);
    
    ShapeF = zeros(nodes_per_element, nqpoints);
    
    for k = 1:nqpoints
        ksi = qPoints(k, 1);
        eta = qPoints(k, 2);
        yip = qPoints(k, 3);
        
        ShapeF(1, k) = 1 - ksi - eta;
        ShapeF(2, k) = ksi;
        ShapeF(3, k) = eta;
        ShapeF(4, k) = yip;
    end
    
elseif size(t, 1) == 7
    
    nodes_per_element = 6;
    
%     Shape function derivative coefficients
    wsp = [1, -3, -3,  2,  4,  2;
        0, -1,  0,  2,  0,  0;
        0,  0, -1,  0,  0,  2;
        0,  4,  0, -4, -4,  0;
        0,  0,  0,  0,  4,  0;
        0,  0,  4,  0, -4, -4];
    
    if numel(nu) == 3 || numel(nu) == nelements*3
    
        qPoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
        weights = [1/6, 1/6, 1/6];
        nqpoints = size(qPoints, 1);
    
    elseif numel(nu) == max(max(t(1:4, :)))*3
        
        qPoints = [1/3, 1/3; 1/5, 1/5; 1/5, 3/5; 3/5, 1/5];
        weights = [-27/96, repmat(25/96, 1, 3)];
        
        viscosity_nodes_per_element = 3;
        nqpoints = size(qPoints, 1);
        ShapeF = zeros(3, nqpoints);
        
        for k = 1:nqpoints
            ksi = qPoints(k, 1);
            eta = qPoints(k, 2);
            
            ShapeF(1, k) = 1 - ksi - eta;
            ShapeF(2, k) = ksi;
            ShapeF(3, k) = eta;
        end
        
    elseif numel(nu) == nnodes*3
        
        qPoints = [0.44594849091597, 0.44594849091597;...
            0.44594849091597, 0.10810301816807;...
            0.10810301816807, 0.44594849091597;...
            0.09157621350977, 0.09157621350977;...
            0.09157621350977, 0.81684757298046;...
            0.81684757298046, 0.09157621350977];
        
        weights = [0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532]/2;
        weights = .5*weights/sum(weights);
        
        viscosity_nodes_per_element = 6;  
        nqpoints = size(qPoints, 1);
        ShapeF = zeros(6, nqpoints);
        
        for i = 1:nodes_per_element
            for k = 1:nqpoints
                ksi = qPoints(k, 1);
                eta = qPoints(k, 2);
                
                coords = [1, ksi, eta, ksi^2, ksi*eta, eta^2]';
                
                ShapeF(i, k) = wsp(i, :)*coords;
            end
        end
    else        
        error('Incorrect diffusivity/viscosity argument')
    end
    
    Nksi = zeros(6, nqpoints);
    Neta = Nksi;
    
    for i = 1:6
        for k = 1:nqpoints
            ksi = qPoints(k, 1);
            eta = qPoints(k, 2);
            
            coords = [zeros(1, 2); eye(2); 2*ksi, 0; eta, ksi; 0, 2*eta];
            
            Nke = wsp(i, :)*coords;
            
            Nksi(i, k) = Nke(1);
            Neta(i, k) = Nke(2);
        end
    end
else
    error('Only 1st and 2nd order meshes are currently supported.')
end

if numel(nu) == 3  % nu is a constant
    const_nu_flag = true;
elseif length(nu(1,:)) == nelements  % nu is defined for each element
    element_nu_flag = true; 
end

indRange = 1:nodes_per_element^2;
KValues = zeros(nodes_per_element^2*nelements, 1);

for el = 1:nelements
    
    J = p(:, t(2:4, el)) - p(:, t(1, el));
    
    detJ = det(J);
    L = inv(J);
    
    Nx = Nksi*L(1, 1) + Neta*L(2, 1) + Nyip*L(3,1);
    Ny = Nksi*L(1, 2) + Neta*L(2, 2) + Nyip*L(3,2);
    Nz = Nksi*L(1, 3) + Neta*L(2, 3) + Nyip*L(3,3);
    
    if const_nu_flag
        K = detJ*(nu(1,:)*weights.*Nx*Nx' + nu(2,:)*weights.*Ny*Ny' + nu(3,:)*weights.*Nz*Nz');
        
    elseif element_nu_flag
        K = detJ*(nu(1,el)*weights.*Nx*Nx' + nu(2,el)*weights.*Ny*Ny' + nu(3,el)*weights.*Nz*Nz');
        
    else
        wnuk = (t(1:viscosity_nodes_per_element, el))*ShapeF.*weights;
        K = detJ*(nu(1,el)*wnuk.*Nx*Nx' + nu(2,el)*wnuk.*Ny*Ny' + nu(3,el)*wnuk.*Nz*Nz');
        
    end
    
    KValues(indRange) = K(:);
    indRange = indRange + nodes_per_element^2;
end

kI = reshape(repmat(t(1:nodes_per_element, :), nodes_per_element, 1), nodes_per_element^2*nelements, 1);
kJ = reshape(repmat(reshape(t(1:nodes_per_element, :), 1, nodes_per_element*nelements), nodes_per_element, 1), nodes_per_element^2*nelements, 1);

% Assemble matrix and right hand side
D = sparse(kI, kJ, KValues, nnodes, nnodes);
F = zeros(nnodes, 1);

end