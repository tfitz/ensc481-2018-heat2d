function [Kcon, Fcon] = assemble_Kcon(ned, nen, nnp, nel_BC, eltype_BC, ...
                          x, y, IEN_BC, ID, quad_rules, hbar, Tinf)
%%
totaldofs = ned*nnp;
Kcon = zeros( totaldofs, totaldofs);
Fcon = zeros( totaldofs, 1);

%% loop over all elements
for e = 1:nel_BC
    
    nen_e = nen(eltype_BC(e));
    nee = ned*nen_e;
    
    % quad rule for each element type
    el_quad = quad_rules{ eltype_BC(e) };
    
    % setup xe ye
    a = 1:nen_e;
    A = IEN_BC(a,e);
    Xe = x(A);
    Ye = y(A);
    
    if eltype_BC(e) == 1
        
        % 2-Node linw
        [Ke, Fe] = el1_ke_con(Xe, Ye, hbar, Tinf, int32(el_quad.nt), el_quad.xi, el_quad.w);
                
    else
        error('Error: unknown element type\n');
        %break
    end
    
    % add to global convective-stiffness matrix
    %a = 1:nen_e;
    %A = IEN_BC(a,e);
    P = ID(1,A);
    Kcon(P,P) = Kcon(P,P) + Ke;
    
    Fcon(P) = Fe;
 
end


