%%% Messing around with Greg's anaerobic growth curve data
%%% Started 10/25/2021

%%% Load in data from spreadsheet
data = xlsread('/Users/peteruhl/OneDrive/Growth-Curves/Anaerobic Growth/CF_otherstrains_anaerobic_GrowthCurves20SEP21.xlsx','Sheet3');

% bacteria name
Efaecalis.dat = data(2:end,1);
Paeruginosa.dat = data(2:end,2);
Calbicans.dat = data(2:end,3);
Sodorifera.dat = data(2:end,4);

tdata = 20:20:48*60;


%%% parameters
r = 0.3554;
Ks = 0.9;
d = 0.076;
k = 1.5;
delta = 0.0002;

p = [r, Ks, d, k, delta];

%%% do optimization here ==================================================
options = optimset('MaxFunEvals',5000,'Display','iter');
% options = optimset('MaxFunEvals',5000);

bdata = Efaecalis.dat;

A = []; b_opt = []; Aeq = []; Beq = [];
lb = zeros(length(p),1);
ub = [30; 5; 3; 1.0; 1];

tic
% [p,fval,flag,output] = fminsearch(@err,p,options,tdata,bdata);
[p,fval,flag,output] = fmincon(@err,p,A,b_opt,Aeq,Beq,lb,ub,[],options,tdata,bdata);
toc



%%% set up and solve ode
tspan = tdata;
y0 = [0.01, 1];

[t, y] = ode15s(@(t,y) rhs(t,y,p), tdata, y0);

J = err(p, tdata, Efaecalis.dat)
p'

% assign solutions
bact = y(:,1);
nut = y(:,2);

%%% Plot stuff ============================================================

hold on; box on
plot(t,bact,"LineWidth",2)
plot(tdata,Efaecalis.dat,'x')
xlabel("Time (minutes)")
ylabel("Optical Density")
legend("Model","Data")


%%% Functions =============================================================

%%% ODE function
function Xp = rhs(t,X,p)

% parameters
r = p(1);
Ks = p(2);
d = p(3);
k = p(4);
delta = p(5);

Xp = zeros(2,1);

x = X(1);
y = X(2);

% ode function
Xp(1) = (r*y)/(Ks + y)*x*(1 - x/k) - d*x;
Xp(2) = -delta*x*y;

end

%%% Error function
function J = err(p,tdata,bdata)
y0 = [0.01, 1];

% solve ode
[t, y] = ode15s(@rhs, tdata, y0, [], p);

% get error
err = bdata - y(:,1);

J = err'*err;

end
