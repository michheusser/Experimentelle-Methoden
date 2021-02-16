clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Modify the number of frames first, using only those frames which are%
% within calibration range                                               %
% 2. Run the Alphasolver over all those pictures only making sure the time
% value has been suitably modified                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Parameters %%%%%%%%%%%%%%%%%
Ti=47;                            %
T_Infinity=25;                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry and material parameters  %
D_hydraulic = 0.02;                 %
lambda_water=0.58;                  %
% Properties of Aluminium           %
Lambda_Al = 238;                    %
Rho_Al = 2700;                      %
cp_Al = 0.88;                       %
a = Lambda_Al/(1000*Rho_Al*cp_Al);  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workpath='D:\Users\student\Desktop\mexBlueFox\'; %folder the m files           %
basename='D:\Users\student\Desktop\mexBlueFox\Groups\'; % folder of pictures  %
groupname='22';                                                                            %
%%%%% Select the frame numbers from where the colorplay starts and it ends
nframe_start=1;                                                                                      %
nframe_end=100;                                                                                      %
cframes=nframe_end-nframe_start+1;                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&

cd(strcat(basename,groupname));
load(groupname);
%%% Setting the time to zero from where it starts
ctvec(1:cframes)=tvec(nframe_start:nframe_end);%-tvec(10); 
cHue_mean(1:cframes)=Hue_mean(nframe_start:nframe_end);
cT_s(1:cframes)=T_s(nframe_start:nframe_end);

cd(basename);
cd(workpath);
for i=1:cframes
    cT_s(i)=polyval(coeff1,cHue_mean(i));  
    % Alpha Solver
    disp('At time in seconds --- '); disp(num2str(ctvec(i)));
    options= optimset('Display','iter');
    Alpha = fminbnd(@xfer,1000,100000,options,ctvec(i),a,cT_s(i),Ti,T_Infinity,Lambda_Al);
    %%%% heat transfer coefficient at every frame computed
    htc(i)=Alpha;
end
%%%%%%%%% Convert to Nusselt Number%%%%%%%
nu=htc*D_hydraulic/lambda_water; 
Nu0=178.9; % Nusselt number on flat surface for this RE and PR numbers
h1=figure(1);
plot(ctvec,nu/Nu0,'-k*');
xlabel('Time [s]');
ylabel('Nu/Nu0 [-]')
title('Nusselt Number ratio, Nu/Nu0');
saveas(h1,[basename groupname '\nusselt_ratio.jpg'])

h2=figure(2);
plot(ctvec,nu,'-b*');
title('Nusselt values [-]');
xlabel('Time [s]');
ylabel('Nu')
saveas(h2,[basename groupname '\nusselt.jpg'])

h3=figure(3);
plot(ctvec,htc,'-b*');
title('Alpha values');
ylabel('alpha [W/m2/K]')
xlabel('Time [s]');
saveas(h3,[basename groupname '\alpha.jpg'])

h4=figure(4);
plot(ctvec,cT_s,'-b*');
title('Surface temperature values');
ylabel('Ts [degC]')
xlabel('Time [s]');
saveas(h4,[basename groupname '\surfaceT.jpg'])
%% write few informations summarizing the test case
fid=fopen([basename groupname '\params_save.txt'],'w');
fprintf(fid,'%s \r\n', 'Summary of the experiment: ');
fprintf(fid,'%s \r\n', ['Date : ' date]);
fprintf(fid,'%s \r\n \r\n',['Group Name : ' groupname]);
fprintf(fid,'%s \r\n', '##########  Material Properties  ##########');
fprintf(fid,'%s \r\n',['Water thermal conductivity : ' num2str(lambda_water) ' W/m/K']);
fprintf(fid,'%s \r\n',['Aluminium - thermal conductivity : ' num2str(Lambda_Al) ' W/m/K']);
fprintf(fid,'%s \r\n',['Aluminium - density : ' num2str(Rho_Al) ' kg/m3']);
fprintf(fid,'%s \r\n \r\n',['Aluminium - specific heat : ' num2str(cp_Al) ' J/kg/K']);
fprintf(fid,'%s \r\n', '########## Operating Conditions  ##########');
fprintf(fid,'%s \r\n',['Initial temperature : ' num2str(Ti) ' degC']);
fprintf(fid,'%s \r\n \r\n',['Fluid temperature : ' num2str(T_Infinity) ' degC']);
fprintf(fid,'%s \r\n', '##########     Reference Data    ##########');
fprintf(fid,'%s \r\n',['Hydraulic Diameter : ' num2str(D_hydraulic) ' m']);
fprintf(fid,'%s \r\n',['Reference Nusselt Number - smooth case : ' num2str(Nu0) ' []']);

fclose(fid);
