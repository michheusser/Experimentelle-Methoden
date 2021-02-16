clc
clear all 
close all
%% 
basename='D:\Users\student\Desktop\mexBlueFox\';
cd(basename);
groupname=char(input('Enter a name for your group here within quotes  '));
mkdir(strcat(basename,'Groups\',groupname))
%%% Take the picture here%%%%%%%%%%%%
%seq=BlueFoxCapture(nframes,exposure_time,gain,triggerMode,whitebalance,xof
%fset,yoffset,width,height);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nframes=100;
exposure_time=40000;
disp('Press any key as soon as the experiment starts....')
pause;
% seq=BlueFoxCapture(nframes,exposure_time,0,0,2,1008,776,172,134);
% seq=BlueFoxCapture(nframes,exposure_time,0,0,2,592,444,216,122);
% seq=BlueFoxCapture(nframes,exposure_time,0,0,2,630,444,90,80); % test crop
% seq=BlueFoxCapture(nframes,exposure_time,0,0,2,1024,776,160,142);

%seq=BlueFoxCapture(nframes,exposure_time,0,0,2,692,770,136,100); % final version
seq=BlueFoxCapture(nframes,exposure_time,0,0,2,424,446,160,80); % final version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ti=input('Enter initial temperature Ti(¡C)=');                  %
% t=input('Enter corresponding time in secondes t(s)=');          %
Ti=46.5;                                                          %
T_Infinity = 25.0;                                                %
D_hydraulic = 0.02;                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Properties of Aluminium                                         %
Lambda_Al = 238;                                                  %
Rho_Al = 2700;                                                    %
cp_Al = 0.88;                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = Lambda_Al/(1000*Rho_Al*cp_Al);
cd(strcat(basename,'Groups\',groupname));
for i=1:nframes
    Tlc_rgb{i}= (seq(1:3,:,:,i));
    Tlc_rgb{i}= permute(Tlc_rgb{i},[2,3,1]);
    
    % Interchange the red and blue values because of some mistake in the camera
    temp(:,:)=Tlc_rgb{i}(:,:,1);
    Tlc_rgb{i}(:,:,1)=Tlc_rgb{i}(:,:,3);
    Tlc_rgb{i}(:,:,3)=temp(:,:);
    
    s=size(Tlc_rgb{i});
    filename=strcat(groupname,num2str(i),'.jpg');
    imwrite(Tlc_rgb{i},filename,'jpg');
    %%% To show picture or not
    %imshow(Tlc_rgb{i},[0,255])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%Time Vector create%%%%%%
t=input('Enter average elapsed time per frame in seconds (s)=');
tvec(1:nframes)=[0:t:t*(nframes-1)];

%% RGB to HSI conversion
for i=1:nframes
    Tlc_hsi{i}=rgb2hsv(double(Tlc_rgb{i}));
    n=1;
    for j=1:s(1)
        for k=1:s(2)
            if mean(Tlc_rgb{i}(j,k,:)) < 250 % making sure there are no extremely high values (in H, S or I) due
                    %to errors in photography, no outliers
                allhue(i,n)=Tlc_hsi{i}(j,k,1);
                n=n+1;
            end
        end
    end
    Hue_mean(i)=mean(allhue(i,:));
end

%% %%%%Calculate the curve fit coeff1 from the calibration pics
cd(basename);
curvefithsi_BL;
%%%%%%%%%%%%IF CURVEFITHSI DOES NOT WORK%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%LOAD THE COEFFICIENTS DIRECTLY%%%%%%%%%%%%%%%%%%
% cd(strcat(basename,'htcalibrate2009\'));
% load('coeff1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% come back to the work directory
cd(strcat(basename,'Groups\',groupname));
    
for i=1:nframes
    T_s(i)=polyval(coeff1,Hue_mean(i));    
    
    % Alpha Solver
    % disp('At time in seconds --- '); disp(num2str(tvec(i)));
    % options= optimset('Display','iter');
    % Alpha =
    % fminbnd(@xfer,1000,100000,options,tvec(i),a,T_s(i),Ti,T_Infinity,Lambda_Al)
    
    %%%% heat transfer coefficient at every frame computed
    % htc(i)=Alpha;  
    %%%%%%%%%
end
save(groupname,'tvec','Hue_mean','T_s','coeff1','nframes');
%%%%%Plot the calibration Curve%%%%%%%%%
h=figure;
plot(e_h,T_cont, '-r',hue_meanc,T, '-k*');
title('Temperature as a Function of Hue-value');
ylabel('Temperature [degC]'); xlabel('Hue-value'); grid on;
legend('interpolated data','raw data');
saveas(h,[basename '\Groups\' groupname '\calibration.jpg'])
% figure(2);
% plot(tvec,htc,'-k*');
% ylabel('Heat Transfer Coefficient');
%%%1. Convert HTC to Nu
%%%2. Discard values outside 
    