%%% Toy model of growth in a chemostat for Greg
%%%
%%% started 9/29/22

global rmax ks dmax ksp lambda mu eta

rmax = 5;
ks = 1;
dmax = 2;
ksp = 0.1;
lambda = 120;
mu = 100;
eta = 1e-8;

B0 = 1e5;
x0 = lambda/mu;

y0 = [B0; x0];
tspan = [0 16];

[t,y] = ode15s(@(t,y) rhs(t,y), tspan, y0);




%%% Plots =================================================================

figure()
hold on; box on;
plot(t,log10(y(:,1)),'linewidth',2)
xlabel('Time','fontsize',16)
legend('Bacteria','fontsize',16)

figure() 
hold on; box on;
plot(t,y(:,2),'linewidth',2)
xlabel('Time','fontsize',16)
legend('Sugar','fontsize',16)


%%% Funcitons =============================================================

function yp = rhs(t,y)
global rmax ks dmax ksp lambda mu eta

B = y(1);
x = y(2);

yp(1) = (rmax*x/(ks + x) - dmax*(1 - x/(ksp + x)))*B;
yp(2) = lambda - mu*x - eta*B*x;
yp = yp';

end