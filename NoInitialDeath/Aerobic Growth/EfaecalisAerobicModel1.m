%%% E faecalis aerobic growth, model #1
%%% model is
%%%
%%% x' = r*x*(1 - (x+y)/k) - dx
%%% y' = dx - gamma*y
%%%
%%% started 2/9/22

% close all

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/Aerobic Growth';
% sheet = sheet + "/Calbicans Aerobic.xlsx";
sheet = sheet + "/Efaecalis Aerobic.xlsx";

%%% Load in data from spreadsheet
data = xlsread(sheet);

% throw out first data point
tdata = data(2:end,1)/60/24;
Ef = data(2:end,2);

%%% parameters
r = 49.9;
d = .3519;
gamma = 20.4378;
alpha_bar = .4;

% IC for ODE
b0 = 0.000384;

p = [r, d, gamma, alpha_bar, b0];

%%% do optimization here ==================================================
options = optimset('MaxFunEvals',10000,'Display','iter'); %,...
%                    'Tolx',1e-1,'TolCon',1e-1,'TolFun',1e-1);


bdata = Ef;

% A = []; b_opt = []; 
Aeq = []; Beq = [];

% this makes sure that gamma > mu
A = [0 0 0 0 0]; b_opt = 0;
lb = zeros(length(p),1);
ub = [50; 50; 50; 100; 0.01];


tic
[p,fval,flag,output] = fminsearch(@err,p,options,tdata,bdata);
% [p,fval,flag,output] = fmincon(@err,p,A,b_opt,Aeq,Beq,lb,ub,[],options,tdata,bdata);
toc

% global tdata bdata
% p = patternsearch(@err,p,A,b,Aeq,Beq,lb,ub,options)


%%% set up and solve ode ==================================================
tspan = tdata;
b0 = p(end);
y0 = [b0, 0];

[t, y] = ode15s(@(t,y) rhs(t,y,p), tdata, y0);

J = err(p, tdata, Ef)
p = p'

%%% calculate AIC
M = length(Ef); % number of data points
Np = length(p); % number of parameters
aic = M*log(J/M) + (2*M*(Np + 1))/(M - Np - 2)

%%% Plot stuff ============================================================
% alpha = 1;
alpha_bar = p(end-1);

figure()
hold on; box on
plot(t,alpha_bar*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Ef,'x')
plot(t,alpha_bar*y(:,1),'Linewidth',2)
plot(t,alpha_bar*y(:,2),'Linewidth',2)
xlabel("Time (days)")
ylabel("Optical Density")
legend("Model","Data","Living","Dead","Location","Northwest")

figure()
hold on; box on
plot(t,alpha_bar*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Ef,'x')
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
d = p(2);
gamma = p(3);

Xp = zeros(2,1);

x = X(1);
y = X(2);

% ode function
Xp(1) = r*x*(1 - (x + y)) - d*x;
Xp(2) = d*x - gamma*y;

end

%%% Error function
function J = err(p,tdata,bdata)
% global tdata bdata

alpha = p(end-1);

b0 = p(end);

y0 = [b0, 0];

% solve ode
[t, y] = ode15s(@rhs, tdata, y0, [], p);

% get error
err = bdata - alpha*((y(:,1) + y(:,2)));

J = err'*err;

end