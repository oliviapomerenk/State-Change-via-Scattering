% evaluate RHS of (22) for a given n at point x(i)
function R = RHS(n,i,n0,k0,A,x,w,N,k,lambda,n_crit)
    integrand = zeros(N,1);
    for j=1:N
        integrand(j) = Gn(n,x(i)-x(j),k,lambda,n_crit) * A(n,n0,x(j)) * exp(1i*k0*x(j));
    end
    R = -sum(integrand .* w);
end
