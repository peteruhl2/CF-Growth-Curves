%%% bar graphs for new anerobic and anaerobic growth for all bugs
%%% 12/21/22

% sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath';
data = xlsread('/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AnaerobicNID/Calbicans Anaerobic.xlsx');

tdata_an = data(2:end,1)/60/24;

%%% second bar graph weight
w2 = 0.5;


%%% C albicans ============================================================
cal_an = data(2:end,2);

% sheet = sheet + "/Calbicans Aerobic.xlsx";
data = xlsread("/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AerobicNID/AlbicansAero.xlsx");
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
axis([0 tdata_an(end) 0 1.5])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18,'location','northwest')
exportgraphics(gcf,'bar_cal.pdf','ContentType','vector')

%%% E faecalis ============================================================
data = xlsread('/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AnaerobicNID/Efaecalis Anaerobic.xlsx');
Ef_an = data(2:end,2);

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AerobicNID/EfaecalisAero.xlsx';
% sheet = sheet + "/Efaecalis Aerobic.xlsx";
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
axis([0 tdata_an(end) 0 1.5])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18,'location','northwest')
exportgraphics(gcf,'bar_ef.pdf','ContentType','vector')

%%% P aeruginosa ============================================================
data = xlsread('/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AnaerobicNID/Paeruginosa Anaerobic.xlsx');
pae_an = data(2:end,2);

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AerobicNID/PseudoAero.xlsx';
% sheet = sheet + "/Paeruginosa Aerobic.xlsx";
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
axis([0 tdata_an(end) 0 1.5])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18,'location','northwest')
exportgraphics(gcf,'bar_pae.pdf','ContentType','vector')

%%% S odorifera ============================================================
data = xlsread('/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AnaerobicNID/Sodorifera Anaerobic.xlsx');
sod_an = data(2:end,2);

sheet = '/Users/peteruhl/Library/CloudStorage/OneDrive-Personal/Growth-Curves/CF-Growth-Curves/NoInitialDeath/AerobicNID/OdorAero.xlsx';
% sheet = sheet + "/Sodorifera Aerobic.xlsx";
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
axis([0 tdata_an(end) 0 1.5])
legend('Anaerobic', 'Aerobic', 'Fontsize', 18,'location','northwest')
exportgraphics(gcf,'bar_sod.pdf','ContentType','vector')