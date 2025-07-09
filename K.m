% calculate kernel of LHS integrals for n1 and n2 at x(i) and x(j)
function ker = K(n1,n2,i,j,A,x,k,lambda,n_crit)
    ker = Gn(n1,x(i)-x(j),k,lambda,n_crit) * A(n1,n2,x(j));
end