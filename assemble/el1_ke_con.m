function [ke, fe] = el1_ke_con( xe, ye, hbar, Tinf, quad_nt, quad_xi, quad_w) %#codegen
% Stiffness matrix, from a convective BC, for a 2-node line
% uniform properties

ned = 1;
nen_e = 2;
nee = ned*nen_e;

%% initialize the matricies
ke = zeros(nee,nee);
fe = zeros(nee,1);

%% Define Integration loop
for i1 = 1:quad_nt
    
    xi = quad_xi(i1);
    
    %% Get Shape Functions, and local derivatives
    [N,Nxi] = el1_ShapeFunctions(xi);
    
    %% Compute Jacobian
    dxdxi = Nxi*xe;
    dydxi = Nxi*ye;
    detJ = sqrt( dxdxi^2 + dydxi^2 );
    
    %% compute dV0
    dV0 = detJ*quad_w(i1);
    
    %% Build ke
    ke = ke - hbar*(N'*N)*dV0;
    
    %% build fe
    fe = fe + hbar*N'*Tinf*dV0;
    
end

