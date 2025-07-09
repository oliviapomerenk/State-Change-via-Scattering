clear;
clc;
close all

set(groot, 'DefaultLineLineWidth', 2);

n0 = 1; % energy state of particle in the box


% both particles are electrons
m1 = 9.109e-31;
m2 = m1;

hbar = 1.055e-34;

L = 1; % size of box for particle 2

U0 = 1e-21;
l0 = 1e-15;
mu0 = U0*l0;


pn_instances = cell(6,1);


k0 = 20;


% can make this go to 100, but takes a while to run. Faster to just take T
% up to 60
for T = 10:10:60
    N = 2*T;
    [Phi,pn_plus,pn_minus,rel_err,x,n_crit,~] = Solve_20(0,L,m1,m2,hbar,mu0,n0,k0,N,T,0);
    pn_instances{T/10} = pn_plus + pn_minus; % vector of size ncrit - 1
end



%%

to_plot_vals = zeros(5,1);
for k=1:5
    to_plot_vals(k) = norm(pn_instances{k+1}-pn_instances{k});
end


figure;
loglog(20:10:60,to_plot_vals(1:5));
ax = gca;
ax.FontSize = 20;
xlabel('truncation order, $T$', 'Interpreter','latex','FontSize',20);
ylabel('successive error in $p_n$','Interpreter','latex','FontSize',20);
set(gca,'TickLength',[0.02, 0.05]);
box on



b=polyfit(log10(20:10:60),log10(to_plot_vals(1:5)),1);





