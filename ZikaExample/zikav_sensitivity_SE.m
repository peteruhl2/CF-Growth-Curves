function zikav_sensitivity_SE
% This function computes local senstivity derivatives and standard erros for 
% parameter estimates for the zika model from Munsur et. al 
%
% INPUT: parameter etimates
% OUPUT: sensitivity graphs, standard errors
% % x1= sh, x2= eh, x3= ih, x4= em, x5=im  % ode variables (Munsur et. al)

% Kidist B-Maxwell 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------
% Add necessary paths
%-----------------------------


addpath('./ode45c')
addpath('./data')
% close all
clear all

%-----------------------------
h= 1e-40;  % Complex step size
%-----------------------------
% Estimated parameters for an Island 
% Yap
bhat=0.9056 ;
bm= 0.0654;%1.09;
ah=1/6.0015;%1/10.57;
gh =1/11.0865;
lm= 1/15;
am= 1/10;
eta= 0.2059;
%-----------------------------
Nh=6892; % Population size
%-----------------------------
% data
data=xlsread('YN.xlsx');
d=data(:,2);
n=16;

%-----------------------------
% add step size to the parameters 
bhatc=bhat+1i*h;
bmc= bm+1i*h;
ahc= ah+ 1i*h;
ghc= gh + 1i*h;
etac = eta + 1i*h;
%-----------------------------
% solve the system with out perturbation 
    function dx= model(t,x)
        dx = [-bhat*x(1)*x(5); bhat*x(1)*x(5)-ah*x(2); ...
              ah*x(2)-gh*x(3);bm*x(3)*(1-x(4)-x(5))-(lm+am)*x(4); ...
              am*x(4)-lm*x(5); eta*ah*x(2)*Nh];
    end
%-----------------------------
% solve the system with perturbed paramater
   function dx1= sens1(t,x)  % bhat
        dx1 = [-bhatc*x(1)*x(5); bhatc*x(1)*x(5)-ah*x(2); ...
              ah*x(2)-gh*x(3);bm*x(3)*(1-x(4)-x(5))-(lm+am)*x(4); ...
              am*x(4)-lm*x(5); eta*ah*x(2)*Nh];
   end
    function dx2= sens2(t,x)  % bm
        dx2 = [-bhat*x(1)*x(5); bhat*x(1)*x(5)-ah*x(2); ...
              ah*x(2)-gh*x(3);bmc*x(3)*(1-x(4)-x(5))-(lm+am)*x(4); ...
              am*x(4)-lm*x(5); eta*ah*x(2)*Nh];
    end
    function dx3= sens3(t,x)  % ah
        dx3 = [-bhat*x(1)*x(5); bhat*x(1)*x(5)-ahc*x(2); ...
              ahc*x(2)-gh*x(3);bm*x(3)*(1-x(4)-x(5))-(lm+am)*x(4); ...
              am*x(4)-lm*x(5); eta*ahc*x(2)*Nh];
    end
    function dx4= sens4(t,x) % gh 
        dx4 = [-bhat*x(1)*x(5); bhat*x(1)*x(5)-ah*x(2); ...
              ah*x(2)-ghc*x(3);bm*x(3)*(1-x(4)-x(5))-(lm+am)*x(4); ...
              am*x(4)-lm*x(5); eta*ah*x(2)*Nh];
    end
    function dx5= sens5(t,x) % eta
        dx5 = [-bhat*x(1)*x(5); bhat*x(1)*x(5)-ah*x(2); ...
              ah*x(2)-gh*x(3);bm*x(3)*(1-x(4)-x(5))-(lm+am)*x(4); ...
              am*x(4)-lm*x(5); etac*ah*x(2)*Nh];
    end
%-----------------------------
% time step
tspan=0:1:n*7; 
%-----------------------------
% initial conditions
IC= [Nh/Nh; 1/Nh; 1/Nh; .005; .005;  Nh/Nh];
%-----------------------------
% solve ODE at e_i , where e_j is the jth parameter
[t,x]= ode45(@model, tspan, IC);
% solution at e_j + ih, where e_j is the jth  parameter  
[t1,s1]= ode45c(@sens1, tspan, IC);
[t2,s2]= ode45c(@sens2, tspan, IC);
[t3,s3]= ode45c(@sens3, tspan, IC);
[t4,s4]= ode45c(@sens4, tspan, IC);
[t5,s5]= ode45c(@sens5, tspan, IC);

% [t1,s1]= ode45(@sens1, tspan, IC);
% [t2,s2]= ode45(@sens2, tspan, IC);
% [t3,s3]= ode45(@sens3, tspan, IC);
% [t4,s4]= ode45(@sens4, tspan, IC);
% [t5,s5]= ode45(@sens5, tspan, IC);
%-----------------------------

%  compute derivatives using the complex-step formula
% dx/de_j ~ imag(x(t,e_j+ih))/h 

ss1=imag(s1(8:7:end,6))/h;  % bh
ss2=imag(s2(8:7:end,6))/h;  % bm
ss3=imag(s3(8:7:end,6))/h;  % ah 
ss4=imag(s4(8:7:end,6))/h;  % gh
ss5=imag(s5(8:7:end,6))/h;  % eta
%-----------------------------
% sensitivity matrix
 M= [ss1, ss2, ss3, ss4, ss5];
 %-----------------------------
%visualize
figure()
plot(M(:,1),'.');
hold on 
plot(M(:,2),'*');
hold on 
plot(M(:,3),'o');
 
hold on 
plot(M(:,4),'>');
hold on 
plot(M(:,5));
legend('dP/d\beta_h','dP/d\beta_m','dP/d\alpha_h','dP/d\gamma_h','dP/d\eta_h','Location','northwest')
axis tight
title('Yap')

figureHandle=gcf;
set(findall(figureHandle,'type','text'),'fontsize',14)
set(gca,'fontsize',14)

%-----------------------------
% Compute standard error
xw= x(8:7:end,6);
J= sum((xw-d).^2);
sigma2= J/(n-2);
MM=inv(M'*M);
dM= diag(MM);
sd= sqrt(sigma2.*dM);
sd=sd';

SD=[sd(1) sd(2) abs((1/ah)-(1/(ah+sd(3)))) abs((1/gh)-(1/(gh+sd(4)))) sd(5)*100];

end


