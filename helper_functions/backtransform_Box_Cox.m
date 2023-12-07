function f_backtr = backtransform_Box_Cox(x_tr, f_tr, x, alpha)
    % x_tr = grid in transformed domain x_tr = g_func(x)
    % f_tr = density in transformed domain
    % x = target grid
    % g_x = transform function g applied to x

    x_tr = x_tr(:);
    f_tr = f_tr(:);
    x = x(:);
    g_x = Box_Cox(x, alpha);


    chain_rule = abs( d_dx_Box_Cox(x, alpha) ); 
    f_backtr = interp1(x_tr, f_tr, g_x, 'linear', 0) .* chain_rule;
end