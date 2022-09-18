%%% S odorifera anaerobic growth, model #3
%%% model is
%%%
%%% x' = ((r*z^n)/(Ks^n + z^n))*x*(1 - x/k) - dx
%%% y' = dx - gamma*y
%%% z' = -delta*x*z + mu*y
%%%
%%% started 2/9/22

% close all

global k
k = 10^9;

%%% Load in data from spreadsheet
data = xlsread('/Users/peteruhl/OneDrive/Growth-Curves/Anaerobic Growth/Sodorifera Anaerobic.xlsx');

% throw out first data point
tdata = data(2:end,1)/60/24;
Sod = data(2:end,2);

%%% parameters
r = 29.9;
Ks = .016835;
n = 1.33;
d = 5.643;
gamma = 4.78;
delta = .0435;
mu = .047;
alpha = 1.5e-9;

% IC for ODE
b0 = 0.0009;
% b0 = Ef(1)/alpha;

p = [r, Ks, n, d, gamma, delta, mu, alpha, b0];
% p = [r, Ks, n, k, d, gamma, delta, mu, b0];

%%% do optimization here ==================================================
options = optimset('MaxFunEvals',5000,'Display','iter'); %,...
%                    'Tolx',1e-1,'TolCon',1e-1,'TolFun',1e-1);


bdata = Sod;

A = []; b_opt = []; 
Aeq = []; Beq = [];


% this makes sure that gamma > mu
% A = [0 0 0 0 -1 0 1 0 0]; b_opt = [0];
lb = zeros(length(p),1);
% lb(end) = 5e4;
ub = [50; 2; 10; 30; 30; 30; 30; 1e-6; 0.1];


tic
% [p,fval,flag,output] = fminsearch(@err,p,options,tdata,bdata);
[p,fval,flag,output] = fmincon(@err,p,A,b_opt,Aeq,Beq,lb,ub,[],options,tdata,bdata);
toc


%%% set up and solve ode ==================================================
tspan = tdata;
b0 = p(end);
y0 = [b0, 0, 1];

[t, y] = ode15s(@(t,y) rhs(t,y,p), tdata, y0);

J = err(p, tdata, Sod)
p = p'

%%% calculate AIC
M = length(Ef); % number of data points
Np = length(p); % number of parameters
aic = M*log(J/M) + (2*M*(Np + 1))/(M - Np - 2)

%%% Plot stuff ============================================================
% alpha = 1;
alpha = p(end-1);

figure()
hold on; box on
plot(t,k*alpha*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Sod,'x')
plot(t,k*alpha*y(:,1),'Linewidth',2)
plot(t,k*alpha*y(:,2),'Linewidth',2)
xlabel("Time (days)")
ylabel("Optical Density")
legend("Model","Data","Living","Dead")

% figure(); hold on
% thing = max(y(:,1)+y(:,2));
% thing2 = max(Pae);
% 
% plot(t,thing2*(y(:,1)+y(:,2))/thing)
% plot(t,Pae)


%%% Functions =============================================================
% p = [r, Ks, n, k, d, gamma, delta, alpha, b0];
%%% ODE function
function Xp = rhs(t,X,p)
global k

% parameters
r = p(1);
Ks = p(2);
n = p(3);
% k = p(4);
d = p(4);
gamma = p(5);
delta = p(6);
mu = p(7);

Xp = zeros(3,1);

x = X(1);
y = X(2);
z = X(3);

% ode function
Xp(1) = ((r*z^n)/(Ks^n + z^n))*x*(1 - (x + y)) - d*x;
Xp(2) = d*x - gamma*y;
Xp(3) = -delta*x*k*z + mu*k*y;

end

%%% Error function
function J = err(p,tdata,bdata)
global k

alpha = p(end-1);
% alpha = 1;
b0 = p(end);
y0 = [b0, 0, 1];

% solve ode
[t, y] = ode15s(@rhs, tdata, y0, [], p);

% get error
err = bdata - alpha*k*((y(:,1) + y(:,2)));

J = err'*err;

end