%%% Sensitivity matrix and standard deviations for C albicans model 1 (anaerobic)
%%% model is
%%%
%%% x' = ((r*z^n)/(Ks^n + z^n))*x*(1 - (x+y)/k) - dx
%%% y' = dx - gamma*y
%%%
%%% started 1/3/23

%%% complex ode45 path
addpath('./ode45c')

% %%% Load in data from spreadsheet
% data = xlsread('/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/Anaerobic Growth/Efaecalis Anaerobic.xlsx');

%%% Load in data from spreadsheet
sheet = pwd;
sheet = sheet + "/Calbicans Anaerobic.xlsx";
data = xlsread(sheet);

% throw out first data point
tdata = data(2:end,1)/60/24;
Cal = data(2:end,2);
n_timepoints = length(Cal);

% %%% best fitting parameters
% r = 49.80;
% % Ks_bar = .065;
% % n = 1.01;
% d = 0.51;
% gamma = 21.34;
% % delta_bar = 49.33;
% % mu_bar = 6.16;
% alpha_bar = 0.27;
% 
% % IC for ODE
% b0 = 0.012;
% 
% % p = [r, Ks_bar, n, d, gamma, delta_bar, mu_bar, alpha_bar, b0];
% p = [r, d, gamma, alpha_bar, b0];

%%% best fitting params from EfaecalisAnaerobicModel1.m
p = [26.9412
   25.0000
    7.0881
    0.4570
    0.0008];

%%% complex step size
h = 1e-40;

%%% add complex step to parameters
p_r = p;
p_r(1) = p_r(1) + 1i*h;

p_d = p;
p_d(2) = p_d(2) + 1i*h;

p_gamma = p;
p_gamma(3) = p_gamma(3) + 1i*h;

p_alpha = p;
p_alpha(4) = p_alpha(4) + 1i*h;

p_b0 = p;
p_b0(5) = p_b0(5) +1i*h;

%%% solve system with perturbed parameters
tspan = tdata;
b0 = p(end);
y0 = [b0, 0];


[t_r, s_r] = ode45c(@(t,y) rhs(t,y,p_r), tdata, [p_r(end), 0, 1]);
% [t_Ks, s_Ks] = ode45c(@(t,y) rhs(t,y,p_Ks), tdata, [p_Ks(end), 0, 1]);
% [t_n, s_n] = ode45c(@(t,y) rhs(t,y,p_n), tdata, [p_n(end), 0, 1]);
[t_d, s_d] = ode45c(@(t,y) rhs(t,y,p_d), tdata, [p_d(end), 0, 1]);
[t_gamma, s_gamma] = ode45c(@(t,y) rhs(t,y,p_gamma), tdata, [p_gamma(end), 0, 1]);
% [t_delta, s_delta] = ode45c(@(t,y) rhs(t,y,p_delta), tdata, [p_delta(end), 0, 1]);
% [t_mu, s_mu] = ode45c(@(t,y) rhs(t,y,p_mu), tdata, [p_mu(end), 0, 1]);
[t_alpha, s_alpha] = ode45c(@(t,y) rhs(t,y,p_alpha), tdata, [p_alpha(end), 0, 1]);
[t_b0, s_b0] = ode45c(@(t,y) rhs(t,y,p_b0), tdata, [p_b0(end), 0, 1]);

[t, y] = ode45c(@(t,y) rhs(t,y,p), tdata, [p(end), 0, 1]);


%%% get derivatives
alpha_bar = p(end-1);

ss_r = imag(alpha_bar*(s_r(:,1) + s_r(:,2)))/h;
% ss_Ks = imag(alpha_bar*(s_Ks(:,1) + s_Ks(:,2)))/h;
% ss_n = imag(alpha_bar*(s_n(:,1) + s_n(:,2)))/h;
ss_d = imag(alpha_bar*(s_d(:,1) + s_d(:,2)))/h;
ss_gamma = imag(alpha_bar*(s_gamma(:,1) + s_gamma(:,2)))/h;
% ss_delta = imag(alpha_bar*(s_delta(:,1) + s_delta(:,2)))/h;
% ss_mu = imag(alpha_bar*(s_mu(:,1) + s_mu(:,2)))/h;
ss_b0 = imag(alpha_bar*(s_b0(:,1) + s_b0(:,2)))/h;

% add complex step here for alpha since its not in the ODE solve
ss_alpha = imag((alpha_bar + (1i*h))*(s_alpha(:,1) + s_alpha(:,2)))/h;

%%% make sensitiviy matrix
% M = [ss_r, ss_Ks, ss_n, ss_d, ss_gamma, ss_delta, ss_mu, ss_alpha, ss_b0];
M = [ss_r, ss_d, ss_gamma, ss_alpha, ss_b0];

J = sum((alpha_bar*((y(:,1)+y(:,2))) - Cal).^2);
% sigma2 = J/(n_timepoints-2);
sigma2 = J/(n_timepoints - length(p));
MM = inv(M'*M);
dM = diag(MM);
sd = sqrt(sigma2.*dM)
% sd = sd'


%%% plot derivative curves
figure()
hold on; box on;
plot(tdata, ss_r, '--')
% plot(tdata, ss_Ks)
% plot(tdata, ss_n, '-')
plot(tdata, ss_d, 'x')
plot(tdata, ss_gamma, '.-')
% plot(tdata, ss_delta, '*')
% plot(tdata, ss_mu, 'o')
plot(tdata, ss_alpha, '.')
% plot(tdata, ss_b0)
legend('dP/dr','dP/dd','dP/d\gamma',...
       'dP/d\alpha','Location','east','Fontsize',12)
xlabel('Time (days)','Fontsize',18)


figure()
hold on; box on
% plot(tdata, ss_Ks)
plot(tdata, ss_b0)
legend('dP/dx_0','Location','northeast','Fontsize',12)
xlabel('Time (days)','Fontsize',18)












%%% Functions =============================================================

%%% ODE function
function Xp = rhs(t,X,p)

% parameters
r = p(1);
% Ks_bar = p(2);
% n = p(3);
d = p(2);
gamma = p(3);
% delta_bar = p(5);
% mu_bar = p(6);

Xp = zeros(3,1);

x = X(1);
y = X(2);
% z = X(3);

% ode function
Xp(1) = r*x*(1 - (x + y)) - d*x;
Xp(2) = d*x - gamma*y;
% Xp(3) = -delta_bar*x*z + mu_bar*y;

end