function output = inverse_Box_Cox(x, alpha)
    if alpha == 0
        output = exp(x);
    else
        output = (x*alpha + 1).^(1/alpha);
        output = output( (x*alpha + 1)  >=0 ); % to avoid raising negative number to power
    end
    output = output(:);
end

