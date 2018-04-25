function [K] = assemble_K(ned, nen, nnp, nel, eltype, ...
                          x, y, IEN, LM, quad_rules)
%%
totaldofs = ned*nnp;
K = zeros(totaldofs,totaldofs);

%% loop over all elements
for e = 1:nel
    
    nen_e = nen(eltype(e));
    nee = ned*nen_e;
    
    % quad rule for each element type
    el_quad = quad_rules{ eltype(e) };
    
    % setup xe ye
    a = 1:nen_e;
    A = IEN(a,e);
    Xe = x(A);
    Ye = y(A);
    
    if eltype(e) == 3
        
        % 4-Node quad
        [Ke] = el3_ke(Xe, Ye, int32(el_quad.nt), el_quad.xi, el_quad.eta, el_quad.w);
                
    else
        error('Error: unknown element type\n');
        %break
    end
    
    % add to global stiffness matrix
    for loop1 = 1:nee
        i = LM(loop1,e);
        for loop2 = 1:nee
            j = LM(loop2,e);
            K(i,j) = K(i,j) + Ke(loop1,loop2);
        end
    end

end


