function hlist = plot_element_solution(axes_handle, IEN, ID, eltype, x, y, z, qn, varargin)


%% Required inputs
p = inputParser;
addRequired(p, 'axes_handle', @simple_axes_check);
addRequired(p, 'IEN', @isnumeric);
addRequired(p, 'ID', @isnumeric);
addRequired(p, 'eltype', @isnumeric);
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'z', @isnumeric);
addRequired(p, 'qn', @isnumeric);


%% Optional Inputs
addParameter(p,               'e_in', 1:size(IEN,2), @isnumeric);
addParameter(p,         'surf_alpha',           0.8, @isnumeric);
addParameter(p, 'show_surfgridlines',         false, @islogical);
addParameter(p,  'show_element_edge',          true, @islogical);
addParameter(p,         'edge_color',           'k' );
addParameter(p,         'edge_width',           1.2, @isnumeric);

% HD surf options per element
addParameter(p,  'e3_pts_per_edge', 3, @isnumeric);
addParameter(p, 'e10_pts_per_edge', 8, @isnumeric);
addParameter(p,        'smoothing', 1, @isnumeric);


%% Parse the inputs
parse(p, axes_handle, IEN, ID, eltype, x, y, z, qn, varargin{:});

e_in           = p.Results.e_in;
alpha0         = p.Results.surf_alpha;
flag_gridlines = p.Results.show_surfgridlines;
flag_edges     = p.Results.show_element_edge;
edge_color     = p.Results.edge_color;
edge_width     = p.Results.edge_width;

e3_pts_per_edge  = p.Results.e3_pts_per_edge;
e10_pts_per_edge = p.Results.e10_pts_per_edge;
flag_smoothing   = p.Results.smoothing;

ned = size(ID,1);

%%
%hnum = 0;
hlist = [];


if flag_smoothing == 0
    opt_FaceColor = 'flat';
else
    opt_FaceColor = 'interp';
end

%%
for e = e_in
    
    if eltype(e) == 3
        % 4-node quad
        nen_e = 4;
        npoints = e3_pts_per_edge;
        A = IEN(1:nen_e,e);
        idx = ID(1:ned,A);
        ue  = qn(idx);
        xe = x(A);
        ye = y(A);
        ze = z(A);
        
        h1 = plot_el_surface(...
            axes_handle, @el3_ShapeFunctions, xe, ye, ze, ue, ...
            npoints, alpha0, flag_gridlines, edge_color, edge_width, opt_FaceColor, flag_edges);
        hlist = [hlist, h1];
        
        
    elseif eltype(e) == 10
        % 9-Node Quad.
        nen_e = 9;
        npoints = e10_pts_per_edge;
        A = IEN(1:nen_e,e);
        idx = ID(1:ned,A);
        ue  = qn(idx);
        xe = x(A);
        ye = y(A);
        ze = z(A);
        
        h1 = plot_el_surface(...
            axes_handle, @el10_ShapeFunctions, xe, ye, ze, ue, ...
            npoints, alpha0, flag_gridlines, edge_color, edge_width, opt_FaceColor, flag_edges);
        
        hlist = [hlist, h1];
        
    else
        fprintf(2',['Error: unknown element type:',...
            ' eltype(%d) = %d\n'],e,eltype(e));
    end
    
end


%%
function h = plot_el_surface(axes_handle, ShapeFunc, xe,ye,ze,ue,npoints,alpha0,flag_gridlines, edge_color, edge_width, opt_FaceColor, flag_edges)

xi0 = linspace(-1,1,npoints);
CData = zeros([npoints, npoints]);

% face
x1p = zeros(npoints,npoints);
y1p = zeros(npoints,npoints);
z1p = zeros(npoints,npoints);
for i1 = 1:npoints
    for i2 = 1:npoints
        NN = ShapeFunc(xi0(i1),xi0(i2));
        x1p(i1,i2) = NN*xe;
        y1p(i1,i2) = NN*ye;
        z1p(i1,i2) = NN*ze;
        CData(i1,i2) = NN*ue;
    end
end

if flag_gridlines == 0
    h = surface(x1p, y1p, z1p, CData,'Parent',axes_handle,'EdgeColor','none', 'FaceColor', opt_FaceColor);
else
    h = surface(x1p, y1p, z1p, CData,'Parent',axes_handle, 'FaceColor', opt_FaceColor);
end

alpha(h, alpha0);

% draw the edge of the element
if flag_edges
    % add edges
    % line 1
    x1p = zeros(1,npoints-1);
    y1p = zeros(1,npoints-1);
    z1p = zeros(1,npoints-1);
    for i1 = 1:npoints-1
        NN = ShapeFunc(xi0(i1),-1);
        x1p(i1) = NN*(xe);
        y1p(i1) = NN*(ye);
        z1p(i1) = NN*(ze);
    end
    % line 2
    x2p = zeros(1,npoints-1);
    y2p = zeros(1,npoints-1);
    z2p = zeros(1,npoints-1);
    for i1 = 1:npoints-1
        NN = ShapeFunc(1,xi0(i1));
        x2p(i1) = NN*(xe);
        y2p(i1) = NN*(ye);
        z2p(i1) = NN*(ze);
    end
    % line 3
    x3p = zeros(1,npoints-1);
    y3p = zeros(1,npoints-1);
    z3p = zeros(1,npoints-1);
    for i1 = 1:npoints-1
        NN = ShapeFunc(-xi0(i1),1);
        x3p(i1) = NN*(xe);
        y3p(i1) = NN*(ye);
        z3p(i1) = NN*(ze);
    end
    % line 4
    x4p = zeros(1,npoints);
    y4p = zeros(1,npoints);
    z4p = zeros(1,npoints);
    for i1 = 1:npoints
        NN = ShapeFunc(-1,-xi0(i1));
        x4p(i1) = NN*(xe);
        y4p(i1) = NN*(ye);
        z4p(i1) = NN*(ze);
    end
    
    h1 = line([x1p,x2p,x3p,x4p], [y1p,y2p,y3p,y4p], ...
        [z1p,z2p,z3p,z4p],...
        'Color',edge_color,...
        'LineWidth',edge_width,...
        'Parent',axes_handle);
    h = [h,h1];
end
