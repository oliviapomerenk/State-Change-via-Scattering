clear;
clc;
close all

set(groot, 'DefaultLineLineWidth', 2);

plot_vs_k = 1; % 1 to plot against k0, 0 to plot against epsilon

do_born = 1; % toggle to also do Born approximation for comparison

compute_error = 0; % toggle to compute and print the max relative error of Phi (returns NaN if 0)

% if this is set to 1, will try to plot all probabilities, not just n s.t.
% an<0 (but these should be the only ones, no bound states)
show_all_probs = 0;

n0 = 1; % energy state of particle in the box

T = 10; % number of terms to include in truncated sum over n prime
N = 2*T; % number of quadrature points


% both particles are electrons
m1 = 9.109e-31;
m2 = m1;

hbar = 1.055e-34;

L = 300e-12; % 300 picometers

mu0 = 1e-30; % strength of interaction (dimension-ful)
g = m1*L*mu0/hbar^2; % strength of interaction (dimensionless)

num_kvals = 500;

k0L = linspace(.1,30,num_kvals);
kvals = k0L/L;
eps_vals = zeros(num_kvals,1);




%% born

if do_born

    probs_born = cell(num_kvals,1);
    probs_plus_born = probs_born;
    probs_minus_born = probs_born;
    
    for i=1:num_kvals
        k0 = kvals(i);
        eps_vals(i) = m1 * mu0 / (k0 * hbar^2);
        a = @(n) m1/m2 * (pi/L)^2 * (n^2 - n0^2) - k0^2; 
        k = @(n) sqrt(-a(n)); % holds for an < 0
        lambda = @(n) sqrt(a(n)); % holds for an > 0
    
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
        
        % vectors of size n_crit-1 for a given k. These should sum to 1
        pn = zeros(n_crit - 1,1);
        p_plus = pn;
        p_minus = pn;
        for n = 1:(n_crit - 1)
            if n0 == n
                pn_plus = 0; % has to be computed separately
                pn_minus = eps_vals(i)^2 * k0/k(n) * ( sin( ( k0+k(n) )*L/2 )*n0*n / ( (k0+k(n)) * L/2 * (k0^2*L^2/pi^2-n^2) ) )^2;
            elseif mod(n0 + n,2) == 0 % n0 + n is even
                pn_plus = eps_vals(i)^2 * k0/k(n) * ( sin((k0-k(n))*L/2)*n0*n/( (k0-k(n)) * L/2 * (k0^2*L^2/pi^2-n^2) ) )^2;
                pn_minus = eps_vals(i)^2 * k0/k(n) * ( sin((k0+k(n))*L/2)*n0*n/( (k0+k(n)) * L/2 * (k0^2*L^2/pi^2-n^2) ) )^2;
            else
                pn_plus = eps_vals(i)^2 * k0/k(n) * ( cos((k0-k(n))*L/2)*n0*n/( (k0-k(n)) * L/2 * (k0^2*L^2/pi^2-n^2) ) )^2;
                pn_minus = eps_vals(i)^2 * k0/k(n) * ( cos((k0+k(n))*L/2)*n0*n/( (k0+k(n)) * L/2 * (k0^2*L^2/pi^2-n^2) ) )^2;
            end
            p_plus(n) = pn_plus;
            p_minus(n) = pn_minus;
        end
    
        if sum(p_plus + p_minus) > 1
            error('epsilon too high');
        end
        p_plus(n0) = 1 - sum(p_plus + p_minus);
        pn = p_plus + p_minus;
    
        probs_born{i} = pn';
        probs_plus_born{i} = p_plus';
        probs_minus_born{i} = p_minus';
    end
    
    
    probs_mat_born = PadAndTurnToMat(probs_born);
    probs_plus_mat_born = PadAndTurnToMat(probs_plus_born);
    probs_minus_mat_born = PadAndTurnToMat(probs_minus_born);
end

%% full non-approximate solution

probs = cell(num_kvals,1);
probs_plus = probs;
probs_minus = probs;


for i=1:num_kvals
    k0 = kvals(i);
    eps_vals(i) = m1 * mu0 / (k0 * hbar^2);
    [Phi,pn_plus,pn_minus,rel_err,x,n_crit,k] = Solve_20(compute_error,L,m1,m2,hbar,mu0,n0,k0,N,T,show_all_probs);
    pn = pn_plus + pn_minus; % vector of size ncrit - 1
    probs{i} = pn';
    probs_plus{i} = pn_plus';
    probs_minus{i} = pn_minus';
end


probs_mat = PadAndTurnToMat(probs);
probs_plus_mat = PadAndTurnToMat(probs_plus);
probs_minus_mat = PadAndTurnToMat(probs_minus);




%% plotting

if plot_vs_k
    to_plot = k0L;
else
    to_plot = eps_vals;
end

figure(1);
hold on
figure(2);
hold on
figure(3);
hold on
for n=1:size(probs_mat,2)
    % if n == n0
    %     figure(1)
    %     plot(to_plot,probs_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}$']);
    %     figure(2)
    %     plot(to_plot,probs_plus_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}^+$']);
    %     figure(3)
    %     plot(to_plot,probs_minus_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}^-$']);
    % else
    if n ~= n0
        figure(1)
        plot(to_plot,probs_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}$']);
        if do_born
            plot(to_plot,probs_mat_born(:,n),'LineWidth',1,'DisplayName',['$p_{',num2str(n),'}$ Born'],'Color','k');
        end
        figure(2)
        plot(to_plot,probs_plus_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}^+$']);
        if do_born
            plot(to_plot,probs_plus_mat_born(:,n),'LineWidth',1,'DisplayName',['$p_{',num2str(n),'}$ Born'],'Color','k');
        end
        figure(3)
        plot(to_plot,probs_minus_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}^-$']);
        if do_born
            plot(to_plot,probs_minus_mat_born(:,n),'LineWidth',1,'DisplayName',['$p_{',num2str(n),'}$ Born'],'Color','k');
        end
    end
end


figure(1)
ax = gca;
ax.FontSize = 20;
ylabel('probability, $p_n$','Interpreter','latex','FontSize',20);
set(gca,'TickLength',[0.02, 0.05]);
box on
%title(['$n_0=$',' ',num2str(n0)],'Interpreter','latex','Fontsize',25);
hl = legend();
set(hl, 'Interpreter','latex');
if ~plot_vs_k
    xlabel('$\epsilon$','Interpreter','latex','FontSize',20);
    set(gca,'XScale','log');
else
    xlabel('dimensionless wavenumber, $Lk_0$','Interpreter','latex','FontSize',20);
end

figure(2)
ax = gca;
ax.FontSize = 20;
ylabel('transmission probability, $p_n^+$','Interpreter','latex','FontSize',20);
set(gca,'TickLength',[0.02, 0.05]);
box on
%title(['$n_0=$',' ',num2str(n0)],'Interpreter','latex','Fontsize',25);
%hl = legend();
%set(hl, 'Interpreter','latex');
if ~plot_vs_k
    xlabel('$\epsilon$','Interpreter','latex','FontSize',20);
    set(gca,'XScale','log');
else
    xlabel('dimensionless wavenumber, $Lk_0$','Interpreter','latex','FontSize',20);
end

figure(3)
ax = gca;
ax.FontSize = 20;
ylabel('reflection probability, $p_n^-$','Interpreter','latex','FontSize',20);
set(gca,'TickLength',[0.02, 0.05]);
box on
%title(['$n_0=$',' ',num2str(n0)],'Interpreter','latex','Fontsize',25);
%hl = legend();
%set(hl, 'Interpreter','latex');
if ~plot_vs_k
    xlabel('$\epsilon$','Interpreter','latex','FontSize',20);
    set(gca,'XScale','log');
else
    xlabel('dimensionless wavenumber, $Lk_0$','Interpreter','latex','FontSize',20);
end










function mat = PadAndTurnToMat(mycell)
M = max(cellfun(@length, mycell));
mycell = cellfun(@(x) [x, zeros(1, M - numel(x))], mycell, 'un', 0);
mat = cell2mat(mycell);
end
