% calculate Phi, probabilities, relative error, discretized domain
function [Phi,pn_plus,pn_minus,rel_err,x,n_crit,k,p_bound] = Solve_20(compute_error,L,m1,m2,hbar,mu0,n0,k0,N,T,show_all_probs)  
    a = @(n) m1/m2 * (pi/L)^2 * (n^2 - n0^2) - k0^2; 
    k = @(n) sqrt(-a(n)); % holds for an < 0
    lambda = @(n) sqrt(a(n)); % holds for an > 0
    
    
    % a1 should be less than 0
    if a(1) >= 0
        error('a1 > 0');
    end

    % an0 should be less than 0
    if a(n0) >= 0
        error('an0 > 0');
    end

    % kn0 should be equal to k0
    if k(n0) ~= k0
        error('kn0 ~= k0');
    end
    
    % find critical n s.t. an changes sign, i.e. a(n_crit) > 0
    n_crit = 1;
    while 1
        if a(n_crit) == 0
            error('a_n=0 exactly, need to deal with this case after all');
        end
        if a(n_crit) > 0
            break;
        end
        n_crit = n_crit + 1;
    end

 
    
    if a(n_crit) < 0 || a(n_crit - 1) > 0
        error('ncrit incorrect');
    end
    
    if T < n_crit
        error('need more terms in truncated sum (must have more than ncrit)');
    end
    
    
    A = @(n,np,x) sin(n*pi*x/L) * sin(np*pi*x/L) * 4*m1*mu0/(hbar^2*L);
    
    
    % set up system

    [x,w] = trap_quad(N,0,L);

    Kmat = zeros(N*T,N*T);

    %wait = waitbar(0,'Starting matrix build...');

    % construct block matrix
    for n1=1:T
        for n2=1:T
            block = zeros(N,N); % one NxN block at the n1,n2 spot
            for i=1:N
                for j=1:N
                    block(i,j) = K(n1,n2,i,j,A,x,k,lambda,n_crit);
                end
            end
            Kmat((n1-1)*N+1:n1*N, (n2-1)*N+1:n2*N) = block; % place block in matrix

        end
        %wait = waitbar(n1/T,wait,'Matrix build in progress...');
    end

    % compute RHS vector
    f = zeros(N*T,1);
    for n=1:T
        fblock = zeros(N,1);
        for i=1:N
            fblock(i) = RHS(n,i,n0,k0,A,x,w,N,k,lambda,n_crit);
        end
        f((n-1)*N+1:n*N) = fblock;
        %wait = waitbar(n/T,wait,'RHS build in progress...');
    end
    
    % solve and shape Phi. Each column corresponds to a given n, each row to a
    % given x(i) on (0,L)

    B = eye(N*T) + Kmat * sparse(diag(repmat(w,T,1)));
    Phi = B\f;
    Phi = reshape(Phi,N,T);


    % compute probabilities
    
    if show_all_probs == 0
        % only calculate these for n s.t. an<0
        pn_plus = zeros(n_crit-1,1);
        pn_minus = zeros(n_crit-1,1);
    
        for n=1:n_crit-1
            if n == n0
                pn_plus(n) = k(n)/k0 * abs(Phi(end,n) + exp(1i*k0*L))^2;
            else
                pn_plus(n) = k(n)/k0 * abs(Phi(end,n))^2;
            end
    
            pn_minus(n) = k(n)/k0 * abs(Phi(1,n))^2;
        end

    else
        % calculate for all
        pn_plus = zeros(T,1);
        pn_minus = zeros(T,1);
        
        for n=1:n_crit-1
            if n == n0
                pn_plus(n) = k(n)/k0 * abs(Phi(end,n) + exp(1i*k0*L))^2;
            else
                pn_plus(n) = k(n)/k0 * abs(Phi(end,n))^2;
            end
            pn_minus(n) = k(n)/k0 * abs(Phi(1,n))^2;
        end
    end
    
    s = sum(pn_plus + pn_minus); % this should be 1
    p_bound = 1 - s; % bound state?

    fprintf('Sum of outgoing probabilities: %f\n',s);

    if compute_error
        rel_err = CalculateMaxError(Phi,N,T,x,w,A,k,lambda,n_crit,n0,k0);
        fprintf('Max relative error: %.5e\n',rel_err);
    else
        rel_err = NaN;
    end

end