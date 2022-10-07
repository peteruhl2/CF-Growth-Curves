%%% S odorifera aerobic growth, model #2
%%% model is
%%%
%%% x' = ((r*z)/(Ks + z))*x*(1 - (x+y)/k) - dx
%%% y' = dx - gamma*y
%%% z' = -delta*x*z + mu*y
%%%
%%% started 2/9/22

% close all

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/Aerobic Growth';
% sheet = sheet + "/Calbicans Aerobic.xlsx";
sheet = sheet + "/Sodorifera Aerobic.xlsx";

%%% Load in data from spreadsheet
data = xlsread(sheet);

% throw out first data point
tdata = data(2:end,1)/60/24;
Sod = data(2:end,2);

%%% parameters
r = 42.9;
Ks_bar = .25;
d = 3.3519;
gamma = 20.4378;
delta_bar = 1.200835;
mu_bar = .0802661;
alpha_bar = .4;

% IC for ODE
b0 = 0.000384;

p = [r, Ks_bar, d, gamma, delta_bar, mu_bar, alpha_bar, b0];

%%% do optimization here ==================================================
options = optimset('MaxFunEvals',10000,'Display','iter'); %,...
%                    'Tolx',1e-1,'TolCon',1e-1,'TolFun',1e-1);


bdata = Sod;

% A = []; b_opt = []; 
Aeq = []; Beq = [];

% this makes sure that gamma > mu
A = [0 0 0 -1 0 1 0 0]; b_opt = 0;
lb = zeros(length(p),1);
ub = [50; 100; 50; 50; 50; 50; 100; 0.01];


tic
% [p,fval,flag,output] = fminsearch(@err,p,options,tdata,bdata);
[p,fval,flag,output] = fmincon(@err,p,A,b_opt,Aeq,Beq,lb,ub,[],options,tdata,bdata);
toc

% global tdata bdata
% p = patternsearch(@err,p,A,b,Aeq,Beq,lb,ub,options)


%%% set up and solve ode ==================================================
tspan = tdata;
b0 = p(end);
y0 = [b0, 0, 1];

[t, y] = ode15s(@(t,y) rhs(t,y,p), tdata, y0);

J = err(p, tdata, Sod)
p = p'

%%% calculate AIC
M = length(Sod); % number of data points
Np = length(p); % number of parameters
aic = M*log(J/M) + (2*M*(Np + 1))/(M - Np - 2)

%%% Plot stuff ============================================================
% alpha = 1;
alpha_bar = p(end-1);

figure()
hold on; box on
plot(t,alpha_bar*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Sod,'x')
plot(t,alpha_bar*y(:,1),'Linewidth',2)
plot(t,alpha_bar*y(:,2),'Linewidth',2)
xlabel("Time (days)")
ylabel("Optical Density")
legend("Model","Data","Living","Dead","Location","Northwest")

figure()
hold on; box on
plot(t,alpha_bar*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Sod,'x')
xlabel("Time (days)", 'Fontsize',18)
ylabel("Optical Density", 'Fontsize',18)
legend("Model","Data",'Fontsize',18,"Location","Northwest")

% figure(); hold on
% thing = max(y(:,1)+y(:,2));
% thing2 = max(Pae);
% 
% plot(t,thing2*(y(:,1)+y(:,2))/thing)
% plot(t,Pae)


%%% Functions =============================================================

%%% ODE function
function Xp = rhs(t,X,p)

% parameters
r = p(1);
Ks_bar = p(2);
d = p(3);
gamma = p(4);
delta_bar = p(5);
mu_bar = p(6);

Xp = zeros(3,1);

x = X(1);
y = X(2);
z = X(3);

% ode function
Xp(1) = ((r*z)/(Ks_bar + z))*x*(1 - (x + y)) - d*x;
Xp(2) = d*x - gamma*y;
Xp(3) = -delta_bar*x*z + mu_bar*y;

end

%%% Error function
function J = err(p,tdata,bdata)
% global tdata bdata

alpha = p(end-1);

b0 = p(end);

y0 = [b0, 0, 1];

% solve ode
[t, y] = ode15s(@rhs, tdata, y0, [], p);

% get error
err = bdata - alpha*((y(:,1) + y(:,2)));

J = err'*err;

end