function [quad_rule] = quad_rules_triangle(opt)

%%
if strcmpi(opt.method, 'prod')
    
    quad_rule = quad_GL_triangle(opt.points);   
    
else
    
    error('Method <%s> not known', opt.method);
    
end
    
    
