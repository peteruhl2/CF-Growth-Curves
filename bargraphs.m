%%% bar graphs for anerobic and anaerobic growth for all bugs
%%% 7/5/22

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/Aerobic Growth';
data = xlsread('/Users/peteruhl/OneDrive/Growth-Curves/Anaerobic Growth/Calbicans Anaerobic.xlsx');

tdata_an = data(2:end,1)/60/24;

%%% second bar graph weight
w2 = 0.5;


%%% C albicans ============================================================
cal_an = data(2:end,2);

sheet = sheet + "/Calbicans Aerobic.xlsx";
data = xlsread(sheet);
cal_aer = data(2:end,2);

tdata_aer = data(2:end,1)/60/24;

% dif = length(cal_an)-length(cal_aer);
% cal_aer = [cal_aer; zeros(dif,1)];

figure()
hold on; box on
bar(tdata_an,cal_an)
bar(tdata_aer,cal_aer, w2)
xlabel("Time (days)", 'Fontsize',18)
ylabel("Optical Density", 'Fontsize',18)
axis([0 tdata_an(end) 0 1.0])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18)
exportgraphics(gcf,'cal_bar.pdf','ContentType','vector')

%%% E faecalis ============================================================
data = xlsread('/Users/peteruhl/OneDrive/Growth-Curves/Anaerobic Growth/Efaecalis Anaerobic.xlsx');
Ef_an = data(2:end,2);

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/Aerobic Growth';
sheet = sheet + "/Efaecalis Aerobic.xlsx";
data = xlsread(sheet);
Ef_aer = data(2:end,2);

tdata_aer = data(2:end,1)/60/24;

figure()
hold on; box on
bar(tdata_an,Ef_an)
bar(tdata_aer,Ef_aer, w2)
% plot(tdata_an,Ef_an,'x')
% plot(tdata_aer,Ef_aer,'x')
xlabel("Time (days)", 'Fontsize',18)
ylabel("Optical Density", 'Fontsize',18)
axis([0 tdata_an(end) 0 1.0])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18)
exportgraphics(gcf,'ef_bar.pdf','ContentType','vector')

%%% P aeruginosa ============================================================
data = xlsread('/Users/peteruhl/OneDrive/Growth-Curves/Anaerobic Growth/Paeruginosa Anaerobic.xlsx');
pae_an = data(2:end,2);

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/Aerobic Growth';
sheet = sheet + "/Paeruginosa Aerobic.xlsx";
data = xlsread(sheet);
pae_aer = data(2:end,2);

tdata_aer = data(2:end,1)/60/24;

figure()
hold on; box on
bar(tdata_an,pae_an)
bar(tdata_aer,pae_aer, w2)
% plot(tdata_an,Ef_an,'x')
% plot(tdata_aer,Ef_aer,'x')
xlabel("Time (days)", 'Fontsize',18)
ylabel("Optical Density", 'Fontsize',18)
axis([0 tdata_an(end) 0 1.0])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18)
exportgraphics(gcf,'pae_bar.pdf','ContentType','vector')

%%% S odorifera ============================================================
data = xlsread('/Users/peteruhl/OneDrive/Growth-Curves/Anaerobic Growth/Sodorifera Anaerobic.xlsx');
sod_an = data(2:end,2);

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/Aerobic Growth';
sheet = sheet + "/Sodorifera Aerobic.xlsx";
data = xlsread(sheet);
sod_aer = data(2:end,2);

tdata_aer = data(2:end,1)/60/24;

figure()
hold on; box on
bar(tdata_an,sod_an)
bar(tdata_aer,sod_aer, w2)
% plot(tdata_an,Ef_an,'x')
% plot(tdata_aer,Ef_aer,'x')
xlabel("Time (days)", 'Fontsize',18)
ylabel("Optical Density", 'Fontsize',18)
axis([0 tdata_an(end) 0 1.0])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18)
exportgraphics(gcf,'sod_bar.pdf','ContentType','vector')
