%%% P aeruginosa anaerobic growth, model #2
%%% model is
%%%
%%% x' = ((r*z)/(Ks + z))*x*(1 - (x+y)/k) - dx
%%% y' = dx - gamma*y
%%% z' = -delta*x*z + mu*y
%%%
%%% started 12/13/22

% close all

global tmax

%%% Load in data from spreadsheet
sheet = pwd;
sheet = sheet + "/Paeruginosa Anaerobic.xlsx";
data = xlsread(sheet);

% throw out first data point
tdata = data(2:end,1)/60/24;
Pae = data(2:end,2);

%%% end of exponential phase
tmax = 0.5;

%%% parameters
r = 22.69;
Ks_bar = .5324475;
d = 1.2925855;
% d = 1.6;
gamma = 3.74378;
delta_bar = 2.509700835;
mu_bar = 7.317502661;
alpha_bar = 0.5852;

% IC for ODE
b0 = 0.00652284;

p = [r, Ks_bar, d, gamma, delta_bar, mu_bar, alpha_bar, b0];

%%% do optimization here ==================================================
options = optimset('MaxFunEvals',10000,'Display','iter'); %,...
%                    'Tolx',1e-1,'TolCon',1e-1,'TolFun',1e-1);


bdata = Pae;

% A = []; b_opt = []; 
Aeq = []; Beq = [];

% this makes sure that gamma > mu
A = [0 0 0 -1 0 1 0 0]; b_opt = 0;
lb = zeros(length(p),1);
ub = [50; 20; 30; 30; 30; 30; 100; 0.01];


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
% [t, y] = ode45c(@(t,y) rhs(t,y,p), tdata, y0);


J = err(p, tdata, Pae)
p = p'

%%% calculate AIC
M = length(Pae); % number of data points
Np = length(p); % number of parameters
aic = M*log(J/M) + (2*M*(Np + 1))/(M - Np - 2)

%%% Plot stuff ============================================================
% alpha = 1;
alpha_bar = p(end-1);

figure()
hold on; box on
plot(t,alpha_bar*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Pae,'x')
plot(t,alpha_bar*y(:,1),'Linewidth',2)
plot(t,alpha_bar*y(:,2),'Linewidth',2)
xlabel("Time (days)")
ylabel("Optical Density")
legend("Model","Data","Living","Dead")

figure()
hold on; box on
plot(t,alpha_bar*(y(:,1)+y(:,2)),"LineWidth",2)
plot(tdata,Pae,'x')
xlabel("Time (days)", 'Fontsize',18)
ylabel("Optical Density", 'Fontsize',18)
legend("Model","Data",'Fontsize',18,'location','east')

figure()
hold on; box on
plot(t,y(:,3), 'Linewidth',2)
xlabel("Time (days)", 'Fontsize',18)
legend("Nutrient",'Fontsize',18,'location','east')


% bar graph of data
% figure()
% hold on; box on
% bar(tdata,Cal)
% xlabel("Time (days)", 'Fontsize',18)
% ylabel("Optical Density", 'Fontsize',18)
% axis([0 tdata(end) 0 1.0])
% legend('C. ablicans optical density', 'Fontsize', 18)
% exportgraphics(gcf,'calan.pdf','ContentType','vector')

% figure(); hold on
% thing = max(y(:,1)+y(:,2));
% thing2 = max(Pae);
% 
% plot(t,thing2*(y(:,1)+y(:,2))/thing)
% plot(t,Pae)


%%% Functions =============================================================

%%% ODE function
function Xp = rhs(t,X,p)
global tmax

% parameters
r = p(1);
Ks_bar = p(2);
% n = p(3);
% d = p(4);
gamma = p(4);
delta_bar = p(5);
mu_bar = p(6);

%%% set no inital death
if t < tmax
    d = 0;
else
    d = p(3);
end

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