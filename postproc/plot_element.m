function hlist = plot_element(axes_handle, e_in, IEN, eltype, x, y, z, varargin)

%% Required inputs
p = inputParser;
addRequired(p, 'axes_handle', @simple_axes_check);
addRequired(p, 'e_in', @isnumeric);
addRequired(p, 'IEN', @isnumeric);
addRequired(p, 'eltype', @isnumeric);
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'z', @isnumeric);

%% Optional Inputs
addParameter(p, 'surf_color', [0.4660,0.6740,0.1880]);
addParameter(p, 'surf_alpha', 0.8, @isnumeric);
addParameter(p, 'edge_color', 'k');
addParameter(p, 'edge_width', 1.2, @isnumeric);


%% Parse the inputs
parse(p, axes_handle, e_in, IEN, eltype, x, y, z, varargin{:});
surf_color = p.Results.surf_color;
alpha0 = p.Results.surf_alpha;
edge_color = p.Results.edge_color;
edge_width = p.Results.edge_width;


%%
% Start the output list
hlist = [];


%%
for e = e_in
    
    if eltype(e) == 2
        %  2: 3-node triangle
        face(1,:) = [1 2 3];
        
        xpts = x(IEN(1:3,e));
        ypts = y(IEN(1:3,e));
        zpts = z(IEN(1:3,e));
        
        h1 = ...
            patch('Parent',axes_handle,'Vertices',[xpts, ypts, zpts],'Faces',face,...
            'FaceColor',surf_color,'FaceAlpha',alpha0, 'EdgeColor', edge_color, 'LineWidth', edge_width);
        hlist = [hlist,h1];
        
        
    elseif eltype(e) == 3
        % 3: 4-node quad
        
        face(1,:) = [1 2 3 4];
        
        xpts = x(IEN(1:4,e));
        ypts = y(IEN(1:4,e));
        zpts = z(IEN(1:4,e));
        
        h1 = ...
            patch('Parent',axes_handle,'Vertices',[xpts, ypts, zpts],'Faces',face,...
            'FaceColor',surf_color,'FaceAlpha',alpha0, 'EdgeColor', edge_color, 'LineWidth', edge_width);
        hlist = [hlist,h1];
        
        
    else
        fprintf(2',...
            'Error: unknown element type: eltype(%d) = %d\n',...
            e,eltype(e));
    end
    
    clear face
end
