clear all;
close all;
clc;

% read temperature & flow
T_1_650=[366, 363, 362];
T_2_650=[272.5, 249.7, 239.6];
T_g_650=[173.2, 176.9, 181.7];
T_1_950=[670.5, 664, 660.5];
T_2_950=[544.9, 490, 463.8];
T_g_950=[292.7, 296.6, 306.6];
flow=[1, 1.5, 2];

% plot temperature
f1 = figure(1)
set(f1,'Position', [0 0 1000 600])
set(f1,'PaperPositionMode','auto')
 
f2 = figure(2)
set(f2,'Position', [0 0 1000 600])
set(f2,'PaperPositionMode','auto')
 
set(0,'CurrentFigure',f1)
s11 = plot(flow,T_1_650,flow,T_2_650,flow,T_g_650)
xlabel('Volumenstrom [nl/min]','FontSize',14)
ylabel('Temperatur [°C]','FontSize',14)
grid on
title('Heiztemperatur 650°C','fontsize', 18)
legend('T_1','T_2','T_{gas}')
 
set(0,'CurrentFigure',f2)
plot(flow,T_1_950,flow,T_2_950,flow,T_g_950)
xlabel('Volumenstrom [nl/min]','FontSize',14)
ylabel('Temperatur [°C]','FontSize',14)
grid on
legend('T_1','T_2','T_{gas}')
title('Heiztemperatur 950°C','fontsize', 18)



%Exporting to LaTeX
print(f1,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/High Temperature/Report/pics/figure1','-dpng')
print(f2,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/High Temperature/Report/pics/figure2','-dpng')

