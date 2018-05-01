%% main
% load GMSH example


%% setup the workspace
clear all
close all
clc

%%
% add in Matlab routines from FEM project
ROOTDIR = fullfile('../..');
path(pathdef)
addpath(fullfile(ROOTDIR,'preprocmesh'));
addpath(fullfile(ROOTDIR,'assemble'));
addpath(fullfile(ROOTDIR,'postproc'));
addpath(fullfile(ROOTDIR,'quadrature'));
addpath(fullfile(ROOTDIR,'shapefunctions'));

%%
% std header stuff
ned = 1;
nen = std_element_defs();

%% Define parameters

% Material constants, Al 2024-T6
alpha = 73e-6; % [m^2/s]
k     = 237  ; % [W/m/K]

% West BC: Dirichlet
T_west = 500;

% East BC: Convection
h    = 30 ; % [ W/m^2/K ]
Tinf = 50 ; % [deg C]
hbar = h/k;

%% Load the msh into matlab
meshfile = 'plate1.msh';
[Nodes, msh, PhysicalName] = load_gmsh(meshfile);

% parse node locations
x = Nodes.x;
y = Nodes.y;
z = Nodes.z;
nnp = length(x);

[nel, eltype, IEN] = parse_msh(msh, 'body');

[nel_BC_s, eltype_BC_s, IEN_BC_s] = parse_msh(msh, 'BC_south');
[nel_BC_e, eltype_BC_e, IEN_BC_e] = parse_msh(msh, 'BC_east');
[nel_BC_n, eltype_BC_n, IEN_BC_n] = parse_msh(msh, 'BC_north');
[nel_BC_w, eltype_BC_w, IEN_BC_w] = parse_msh(msh, 'BC_west');

%% Start plot
[fig,ax] = init_plot_figure();
ax.XLim = [-0.1, 3.1];
ax.YLim = [-0.1, 1.1];

% plot_node(ax, 1:nnp, x,y,z, 'r');
plot_node_labels(ax, 1:nnp, x,y,z);

plot_element(ax, 1:nel, IEN, eltype, x,y,z, 'surf_alpha', 0.2);
plot_element_labels(ax, IEN, eltype, x, y, z);

%% 
% Plot the BC nodes
plot_node(ax, unique(IEN_BC_s(:)), x,y,z, 'r');
plot_node(ax, unique(IEN_BC_e(:)), x,y,z, 'g');
plot_node(ax, unique(IEN_BC_n(:)), x,y,z, 'b');
plot_node(ax, unique(IEN_BC_w(:)), x,y,z, 'm');


%% Apply BC's
% generate BC tables for all nodes

% this is flagged if a node has an essential BC
fix = zeros(ned,nnp);
% this is list of the values of that BC
g_list = zeros(ned,nnp);


% east
% A = unique(IEN_BC_e(:));
% fix(1, A ) = 1;
% g_list(1, A) = T_east;

% west
A = unique(IEN_BC_w(:));
fix(1, A ) = 1;
g_list(1, A) = T_west;


%% generate or load gaussian quadrature information
quad_rules = set_integration_rules( [eltype; eltype_BC_e] );

%% Construct the IM and LM Matricies
[ID, LM, neq, gg, nee, ng, idx_ff, idx_fr] = build_mesh(...
    nnp, IEN, g_list, nen, ned, nel, eltype, fix);

%% build Mass and Stiffness Matrices
ndofs = neq+ng;
qn = zeros(ndofs,1);

%%
% Build M
% [M] = assemble_M(ned, nen, nnp, nel, eltype, ...
%     x, y, IEN, LM, quad_rules, alpha);

%%
% Build K
[K] = assemble_K(ned, nen, nnp, nel, eltype, ...
    x, y, IEN, LM, quad_rules);


%%
% Build Kcon
% convecting on the east side:
[Kcon, Fcon] = assemble_Kcon(ned, nen, nnp, nel_BC_e, eltype_BC_e, ...
                          x, y, IEN_BC_e, ID, quad_rules, hbar, Tinf);
K = K-Kcon;

%% Apply Dirichlet BC's
totaldofs = ned*nnp;
q = zeros(totaldofs, 1);

idx = find(fix==1);
for i = 1:ng % ng = length(idx)
    
    A = idx(i);
    P = ID(1,A);
    
    q(P) = g_list(A);
    
end

%% Solve for the unknowns
K_ff = K(idx_ff,idx_ff);
K_fr = K(1:neq, idx_fr);

R = zeros(neq,1) + Fcon(1:neq);

% solve:
q(1:neq) = K_ff\( R - K_fr*q(idx_fr) );

%% make a new plot to show the solution
[~,ax2] = init_plot_figure();
ax2.XLim = [-0.1, 3.1];
ax2.YLim = [-0.1, 1.1];

plot_element_solution(ax2, IEN, ID, eltype, x, y, z, q);
