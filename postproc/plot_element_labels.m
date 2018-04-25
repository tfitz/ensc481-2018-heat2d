function hlist = plot_element_labels(axes_handle, IEN, eltype, x, y, z,...
    varargin)

%% Parse input
p = inputParser;
addRequired( p, 'axes_handle', @simple_axes_check);
addRequired( p,         'IEN', @isnumeric);
addRequired( p,      'eltype', @isnumeric);
addRequired( p,           'x', @isnumeric);
addRequired( p,           'y', @isnumeric);
addRequired( p,           'z', @isnumeric);

addParameter(p,        'e_in', 1:size(IEN,2), @isnumeric);

addParameter(p,       'field',      'number',  @ischar);
addParameter(p,    'fontsize',            10, @isnumeric);
addParameter(p,   'fontcolor',           'k');
addParameter(p,   'fontweight',       'bold', @ischar);
addParameter(p,'backgroundcolor',     'none');

% Parse the inputs
parse(p, axes_handle, IEN, eltype, x, y, z, varargin{:});

e_in           = p.Results.e_in;

flag_field     = p.Results.field;
fontsize       = p.Results.fontsize;
fontcolor      = p.Results.fontcolor;
fontweight     = p.Results.fontweight;
backgroundcolor= p.Results.backgroundcolor;


%% work by element, placing the label at the center
hlist = nan(size(e_in));

i = 0;
for e = e_in
    i = i+1;
    
    if eltype(e) == 3
        % 4-node quad
        nen_e = 4;
        NN = el3_ShapeFunctions(0,0);
            
    else
        error('Element type <%d> is not yet implemented', eltype(e) );
    end
    
    [xe,ye,ze] = get_xyz(e, nen_e, IEN, x, y, z);
    
    x1p = NN*xe;
    y1p = NN*ye;
    z1p = NN*ze;
    
    if( strcmpi(flag_field, 'number') )
        el_str = sprintf('%d', e);
        
    elseif( strcmpi(flag_field, 'eltype') )
        el_str = sprintf('%d', eltype(e));
        
    else
        error('label type <%s> not implemented', flag_field);
    end
    
    hlist(i) = text(x1p, y1p, z1p, el_str, ...
        'FontSize', fontsize, ...
        'Color', fontcolor, ...
        'BackgroundColor', backgroundcolor, ...
        'FontWeight', fontweight, ...
        'HorizontalAlignment', 'center', ...
        'Parent', axes_handle);
    
end

%%
function [xe, ye, ze] = get_xyz(e, nen_e, IEN, x, y, z)
a = 1:nen_e;
A = IEN(a,e);
xe = x(A);
ye = y(A);
ze = z(A);
