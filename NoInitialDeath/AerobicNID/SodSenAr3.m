%%% Sensitivity matrix and standard deviations for S odorifera model 3
%%% (aerobic)
%%% model is
%%%
%%% x' = ((r*z^n)/(Ks^n + z^n))*x*(1 - (x+y)/k) - dx
%%% y' = dx - gamma*y
%%% z' = -delta*x*z + mu*y
%%%
%%% started 12/31/22

%%% complex ode45 path
addpath('./ode45c')

% %%% Load in data from spreadsheet
% data = xlsread('/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/Anaerobic Growth/Efaecalis Anaerobic.xlsx');

%%% Load in data from spreadsheet
sheet = pwd;
sheet = sheet + "/OdorAero.xlsx";
data = xlsread(sheet);

% throw out first data point
tdata = data(2:end,1)/60/24;
Sod = data(2:end,2);
n_timepoints = length(Sod);

% %%% best fitting parameters
% r = 49.51;
% Ks_bar = .020;
% n = 1.01;
% d = 0.66;
% gamma = 12.04;
% delta_bar = 12.21;
% mu_bar = 0.12;
% alpha_bar = 0.33;
% 
% % IC for ODE
% b0 = 0.0082;
% 
% p = [r, Ks_bar, n, d, gamma, delta_bar, mu_bar, alpha_bar, b0];

%%% best fitting params from EfaecalisAnaerobicModel3.m
p = [52.7646
    1.4962
    1.0298
   10.1645
   27.4253
    7.6789
   25.9565
    2.3549
    0.0050];

%%% complex step size
h = 1e-40;

%%% add complex step to parameters
p_r = p;
p_r(1) = p_r(1) + 1i*h;

p_Ks = p;
p_Ks(2) = p_Ks(2) + 1i*h;

p_n = p;
p_n(3) = p_n(3) + 1i*h;

p_d = p;
p_d(4) = p_d(4) + 1i*h;

p_gamma = p;
p_gamma(5) = p_gamma(5) + 1i*h;

p_delta = p;
p_delta(6) = p_delta(6) + 1i*h;

p_mu = p;
p_mu(7) = p_mu(7) + 1i*h;

p_alpha = p;
p_alpha(8) = p_alpha(8) + 1i*h;

p_b0 = p;
p_b0(9) = p_b0(9) +1i*h;

%%% solve system with perturbed parameters
tspan = tdata;
b0 = p(end);
y0 = [b0, 0, 1];


[t_r, s_r] = ode45c(@(t,y) rhs(t,y,p_r), tdata, [p_r(end), 0, 1]);
[t_Ks, s_Ks] = ode45c(@(t,y) rhs(t,y,p_Ks), tdata, [p_Ks(end), 0, 1]);
[t_n, s_n] = ode45c(@(t,y) rhs(t,y,p_n), tdata, [p_n(end), 0, 1]);
[t_d, s_d] = ode45c(@(t,y) rhs(t,y,p_d), tdata, [p_d(end), 0, 1]);
[t_gamma, s_gamma] = ode45c(@(t,y) rhs(t,y,p_gamma), tdata, [p_gamma(end), 0, 1]);
[t_delta, s_delta] = ode45c(@(t,y) rhs(t,y,p_delta), tdata, [p_delta(end), 0, 1]);
[t_mu, s_mu] = ode45c(@(t,y) rhs(t,y,p_mu), tdata, [p_mu(end), 0, 1]);
[t_alpha, s_alpha] = ode45c(@(t,y) rhs(t,y,p_alpha), tdata, [p_alpha(end), 0, 1]);
[t_b0, s_b0] = ode45c(@(t,y) rhs(t,y,p_b0), tdata, [p_b0(end), 0, 1]);

[t, y] = ode45c(@(t,y) rhs(t,y,p), tdata, [p(end), 0, 1]);


%%% get derivatives
alpha_bar = p(end-1);

ss_r = imag(alpha_bar*(s_r(:,1) + s_r(:,2)))/h;
ss_Ks = imag(alpha_bar*(s_Ks(:,1) + s_Ks(:,2)))/h;
ss_n = imag(alpha_bar*(s_n(:,1) + s_n(:,2)))/h;
ss_d = imag(alpha_bar*(s_d(:,1) + s_d(:,2)))/h;
ss_gamma = imag(alpha_bar*(s_gamma(:,1) + s_gamma(:,2)))/h;
ss_delta = imag(alpha_bar*(s_delta(:,1) + s_delta(:,2)))/h;
ss_mu = imag(alpha_bar*(s_mu(:,1) + s_mu(:,2)))/h;
ss_b0 = imag(alpha_bar*(s_b0(:,1) + s_b0(:,2)))/h;

% add complex step here for alpha since its not in the ODE solve
ss_alpha = imag((alpha_bar + (1i*h))*(s_alpha(:,1) + s_alpha(:,2)))/h;


%%% make sensitiviy matrix
M = [ss_r, ss_Ks, ss_n, ss_d, ss_gamma, ss_delta, ss_mu, ss_alpha, ss_b0];
% M = [ss_r, ss_Ks, ss_n, ss_d, ss_gamma, ss_delta, ss_mu, ss_b0];

J = sum((alpha_bar*((y(:,1)+y(:,2))) - Sod).^2);
% sigma2 = J/(n_timepoints-2);
sigma2 = J/(n_timepoints - length(p));
MM = inv(M'*M);
dM = diag(MM);
sd = sqrt(sigma2.*dM)
% sd = sd'


%%% plot derivative curves
figure()
hold on; box on;
plot(tdata, ss_r, '>')
% plot(tdata, ss_Ks)
plot(tdata, ss_n, '-')
plot(tdata, ss_d, 'x')
plot(tdata, ss_gamma, '.-')
plot(tdata, ss_delta, '*')
plot(tdata, ss_mu, 'o')
plot(tdata, ss_alpha, '.')
% plot(tdata, ss_b0)
legend('dP/dr','dP/dn','dP/dd','dP/d\gamma','dP/d\delta','dP/d\mu',...
       'dP/d\alpha','Location','northwest','Fontsize',12)
xlabel('Time (days)','Fontsize',18)

set(legend,...
    'Position',[0.7535 0.3452 0.1348 0.2583],...
    'FontSize',12);


figure()
hold on; box on
plot(tdata, ss_Ks, '<')
plot(tdata, ss_b0/10 ,'--')
legend('dP/dK_s','dP/dx_0 x 1/10','Location','northeast','Fontsize',12)
xlabel('Time (days)','Fontsize',18)












%%% Functions =============================================================

%%% ODE function
function Xp = rhs(t,X,p)

% parameters
r = p(1);
Ks_bar = p(2);
n = p(3);
d = p(4);
gamma = p(5);
delta_bar = p(6);
mu_bar = p(7);

Xp = zeros(3,1);

x = X(1);
y = X(2);
z = X(3);

% ode function
Xp(1) = ((r*z^n)/(Ks_bar^n + z^n))*x*(1 - (x + y)) - d*x;
Xp(2) = d*x - gamma*y;
Xp(3) = -delta_bar*x*z + mu_bar*y;

end