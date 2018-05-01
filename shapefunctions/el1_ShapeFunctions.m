function [NN,Nr] = el1_ShapeFunctions(r)

%% Compute shape functions for 2 node line

%%
NN = [(1-r)/2, (1+r)/2];

Nr = [-1/2, 1/2];

