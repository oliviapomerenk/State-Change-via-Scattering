clear;
clc;
close all

set(groot, 'DefaultLineLineWidth', 2);

plot_vs_k = 1; % 1 to plot against k0, 0 to plot against epsilon

compute_error = 0; % toggle to compute and print the max relative error of Phi (returns NaN if 0)

% if this is set to 1, will try to plot all probabilities, not just n s.t.
% an<0 (but these should be the only ones, no bound states)
show_all_probs = 0;

n0 = 1; % energy state of particle in the box

T = 50; % number of terms to include in truncated sum over n prime
N = 2*T; % number of quadrature points


% both particles are electrons
m1 = 9.109e-31;
m2 = m1;

hbar = 1.055e-34;

L = 300e-12; % 300 picometers

mu0 = 1e-22; % strength of interaction (dimension-ful)
g = m1*L*mu0/hbar^2; % strength of interaction (dimensionless)

num_kvals = 500;

k0L = linspace(.1,30,num_kvals);
kvals = k0L/L;
eps_vals = zeros(num_kvals,1);


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
        figure(2)
        plot(to_plot,probs_plus_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}^+$']);
        figure(3)
        plot(to_plot,probs_minus_mat(:,n),'LineWidth',2,'DisplayName',['$p_{',num2str(n),'}^-$']);
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
