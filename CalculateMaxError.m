% maximum relative error of Phi based on subbing into (22)
function rel_err = CalculateMaxError(Phi,N,T,x,w,A,k,lambda,n_crit,n0,k0)
    rel_err = 0;
    for i=1:N
        for n=1:T
            LHS_22 = Phi(i,n);
            for np=1:T % sum over n prime
                % add up each integral
                int = 0;
                for j=1:N
                    % int = int + Gn(n,x(i)-x(j),k,lambda,n_crit) * A(n,np,x(j)) * Phi(j,np) * w(j);
                    int = int + K(n,np,i,j,A,x,k,lambda,n_crit) * Phi(j,np) * w(j);
                end
                LHS_22 = LHS_22 + int;
            end
            RHS_22 = RHS(n,i,n0,k0,A,x,w,N,k,lambda,n_crit);
            rel_err = max(rel_err, abs((LHS_22-RHS_22)/RHS_22));
        end
    end
end
