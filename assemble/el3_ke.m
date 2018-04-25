function [ke] = el3_ke( xe, ye, quad_nt, quad_xi, quad_eta, quad_w) %#codegen
% Stiffness matrix of 4-node quad
% uniform properties

ned = 1;
nen_e = 4;
nee = ned*nen_e;

%% initialize the matricies
ke = zeros(nee,nee);

%% Define Integration loop
for i1 = 1:quad_nt
    xi = quad_xi(i1);
    eta = quad_eta(i1);
    
    %% Get Shape Functions, and local derivatives
    [~,Nxi,Neta] = el3_ShapeFunctions(xi,eta);
    
    %% Compute Jacobian
    J = [Nxi*xe, Neta*xe;
        Nxi*ye, Neta*ye];
    detJ = det(J);
    
    %% Derivative matrix B
    B = J\[Nxi; Neta];
    % Nx = B(1,:);
    % Ny = B(2,:);
    
    %% compute dV0
    dV0 = detJ*quad_w(i1);
    
    %% Build ke
    ke = ke + B'*B*dV0;
    
end

