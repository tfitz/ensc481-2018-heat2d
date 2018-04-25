function [me] = el3_me( xe, ye, alpha, quad_nt, quad_xi, quad_eta, quad_w) %#codegen
% Mass matrix Element

%%
ned = 1;
nen_e = 4;
nee = 4;

%% initialize the matricies
me = zeros(ned*nen_e,ned*nen_e);

%% Define Integration loop
for i1 = 1:quad_nt
    xi = quad_xi(i1);
    eta = quad_eta(i1);
    
    %% Get Shape Functions, and local derivatives
    [N,Nxi,Neta] = el3_ShapeFunctions(xi,eta);
        
    %% Compute Jacobian
    J = [Nxi*xe, Neta*xe;
         Nxi*ye, Neta*ye];
    detJ = det(J);
    
    %% compute dV0
    dV0 = detJ*quad_w(i1);
    
    %% Build me
    me = me + 1/alpha*(N'*N)*dV0;
    
end
