%% main_01
% load simple example


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

% Material constants
alpha = 1;


%% Make a new FEM mesh
nel = 9;
nnp = 16;
eltype(1:nel) = 3;

alpha(1:nel) = alpha;

% connectivity
%  IEN(a,e) = A;
%  A: global node number
%  a: local node number
%  e: element number
% I'll put it in by element, and transpose:
IEN = [ 1,2,6,5 ;
    2,3,7,6 ;
    3,4,8,7 ;
    5,6,10,9;
    6,7,11,10;
    8,12,11,7;
    9,10,14,13;
    10,11,15,14;
    11,12,16,15]';
x = zeros(nnp, 1);
y = zeros(nnp, 1);
z = zeros(nnp, 1);

x([1,5,9,13]) = 0;
x([2,6,10,14]) = 1/3;
x([3,7,11,15]) = 2/3;
x([4, 8, 12, 16]) = 1;

y(1:4) = 0;
y(5:8) = 1/3;
y(9:12) = 2/3;
y(13:16) = 1;

%% Start plot
[fig,ax] = init_plot_figure();
ax.XLim = [-0.1, 1.1];
ax.YLim = [-0.1, 1.1];

% plot_node(ax, 1:nnp, x,y,z, 'r');
plot_node_labels(ax, 1:nnp, x,y,z);

plot_element(ax, 1:nel, IEN, eltype, x,y,z, 'surf_alpha', 0.2);
plot_element_labels(ax, IEN, eltype, x, y, z);


%% generate or load gaussian quadrature information
quad_rules = set_integration_rules( eltype );


%% Apply BC's
% generate BC tables for all nodes

% this is flagged if a node has an essential BC
fix = zeros(ned,nnp);
fix(1,[1,2,3,4, 8,12,16, 15, 14, 13, 9, 5]) = 1;

% this is list of the values of that BC
g_list = zeros(ned,nnp);
T_south = 0;
T_east  = 50;
T_north = 100;
T_west = 75;
g_list(1,[2,3]) = T_south;
g_list(1,[8,12]) = T_east;
g_list(1,[15,14]) = T_north;
g_list(1,[9,5]) = T_west;
% I'll just take the temp average for the corner nodes:
g_list(1,[1]) = (T_south+T_west)/2;
g_list(1,[4]) = (T_south+T_east)/2;
g_list(1,[16]) = (T_north+T_east)/2;
g_list(1,[13]) = (T_north+T_west)/2;

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

%% Apply BC's
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

R = zeros(neq,1);

% solve:
q(1:neq) = K_ff\( R - K_fr*q(idx_fr) );

%% make a new plot to show the solution
[~,ax2] = init_plot_figure();
ax2.XLim = [-0.1, 1.1];
ax2.YLim = [-0.1, 1.1];

plot_element_solution(ax2, IEN, ID, eltype, x, y, z, q);
