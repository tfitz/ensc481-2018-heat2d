function [M] = assemble_M(ned, nen, nnp, nel, eltype, ...
                          x, y, IEN, LM, quad_rules, alpha)


%%
totaldofs = ned*nnp;

%% storage
% could be a sparse array if we are crafty and concerned about larger problems
M = zeros(totaldofs, totaldofs);

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
        [Me] = el3_me(Xe, Ye, alpha(e), int32(el_quad.nt), el_quad.xi, el_quad.eta, el_quad.w);
                
    else
        fprintf(2','Error: unknown element type\n');
        %break
    end

    %% insert into Global stiffness matrix
    for loop1 = 1:nee
        i = LM(loop1, e);
        for loop2 = 1:nee
            j = LM(loop2,e);
            M(i,j) = M(i,j) + Me(loop1,loop2);
        end
    end

end

