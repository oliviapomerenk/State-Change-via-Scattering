% evaluate Green's function for given n at x
function G = Gn(n,x,k,lambda,n_crit)
    if x ~= 0
        if n < n_crit
            if x < 0
                G = -1/(2*1i*k(n)) * exp(-1i*k(n)*x);
            else
                G = -1/(2*1i*k(n)) * exp(1i*k(n)*x);
            end
        else
            if x < 0
                G = 1/(2*lambda(n)) * exp(lambda(n)*x);
            else
                G = 1/(2*lambda(n)) * exp(-lambda(n)*x);
            end
        end
    else
        if n < n_crit
            G = -1/(2*1i*k(n));
        else
            G = 1/(2*lambda(n)); % not sure if this is ideal
        end
    end

end