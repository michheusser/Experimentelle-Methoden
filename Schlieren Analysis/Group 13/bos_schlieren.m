close all
clear all
clc

%% 
%Wedge Calibration


% Einlesen der Bilder
r=imread('Referenz.tif');
 
a=imread('75.tif');
b=imread('125.tif');
c=imread('175.tif');
d=imread('225.tif');
e=imread('275.tif');
f=imread('325.tif');
g=imread('375.tif');
h=imread('425.tif');
i=imread('475.tif');
j=imread('525.tif');
k=imread('575.tif');
l=imread('625.tif');
m=imread('675.tif');
n=imread('725.tif');
o=imread('775.tif');
p=imread('825.tif');
 
% Wählen des Bildausschnittes
rect = [551, 377, 157, 142];
 
rs=imcrop(r,rect);
 
as=imcrop(a,rect);
bs=imcrop(b,rect);
cs=imcrop(c,rect);
ds=imcrop(d,rect);
es=imcrop(e,rect);
fs=imcrop(f,rect);
gs=imcrop(g,rect);
hs=imcrop(h,rect);
is=imcrop(i,rect);
js=imcrop(j,rect);
ks=imcrop(k,rect);
ls=imcrop(l,rect);
ms=imcrop(m,rect);
ns=imcrop(n,rect);
os=imcrop(o,rect);
ps=imcrop(p,rect);
 
 
% Korrelation
rc=normxcorr2(rs,r);
 
ac=normxcorr2(as,r);
bc=normxcorr2(bs,r);
cc=normxcorr2(cs,r);
dc=normxcorr2(ds,r);
ec=normxcorr2(es,r);
fc=normxcorr2(fs,r);
gc=normxcorr2(gs,r);
hc=normxcorr2(hs,r);
ic=normxcorr2(is,r);
jc=normxcorr2(js,r);
kc=normxcorr2(ks,r);
lc=normxcorr2(ls,r);
mc=normxcorr2(ms,r);
nc=normxcorr2(ns,r);
oc=normxcorr2(os,r);
pc=normxcorr2(ps,r);

%%


%Font sizes for figures
titlesize= 20;
axissize=18;
set(0,'DefaultAxesFontSize',axissize)

%Creation of figure of proper size
f0 = figure(1)
set(f0,'Position', [0 1000 1000 400])
set(f0,'PaperPositionMode','auto')% print() exports the same size as figures are shown
s01=subplot(1,2,1)
set(s01,'Position',[0.005 0.13 0.45 0.78])

%Compound plot
imshow(imcrop(rc+ac+bc+cc+dc+ec+fc+gc+hc+ic+jc+kc+lc+mc+nc+oc+pc,[590, 400, 157, 142]))
title('Compound Visual Display', 'Fontsize',titlesize)

% Koordinaten der Korrelation
[maxval(1),maxind(1)]=max(rc(:));
[maxval(2),maxind(2)]=max(ac(:));
[maxval(3),maxind(3)]=max(bc(:));
[maxval(4),maxind(4)]=max(cc(:));
[maxval(5),maxind(5)]=max(dc(:));
[maxval(6),maxind(6)]=max(ec(:));
[maxval(7),maxind(7)]=max(fc(:));
[maxval(8),maxind(8)]=max(gc(:));
[maxval(9),maxind(9)]=max(hc(:));
[maxval(10),maxind(10)]=max(ic(:));
[maxval(11),maxind(11)]=max(jc(:));
[maxval(12),maxind(12)]=max(kc(:));
[maxval(13),maxind(13)]=max(lc(:));
[maxval(14),maxind(14)]=max(mc(:));
[maxval(15),maxind(15)]=max(nc(:));
[maxval(16),maxind(16)]=max(oc(:));
[maxval(17),maxind(17)]=max(pc(:));
 
 
% Pixelkoordinaten
size_cr=size(ac); % Grösse des cropped image, gleich für alle Bilder
for ii=1:17
    [rows(ii),cols(ii)]=ind2sub(size_cr,maxind(ii));
end
 
 
% Verschiebung in Pixelkoordinaten
for ii=1:16
    delta(ii,:)=norm([rows(ii+1)-rows(1),cols(ii+1)-cols(1)]);
end
 
 
% Calibration plot creation and exporting
s02 = subplot(1,2,2)
set(s02,'Position',[0.512  0.13 0.45 0.78])
distance=[75,125,175,225,275,325,375,425,475,525,575,625,675,725,775,825];
plot(distance,delta,'-s')
title('Wedge Calibration', 'Fontsize',titlesize)
xlabel('Object Distance [mm]', 'Fontsize',axissize)
ylabel('Apparent Displacement [pix]', 'Fontsize',axissize)

grid on
print(f0,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/Schlieren Analysis/Report/pics/figure0','-dpng')

%% 

% Flow cell

close all
clc

%Creation of figures of the proper size

f1 = figure(2); %Zeitpunkt 1
set(0,'DefaultAxesFontSize',axissize)
set(f1,'Position',[0 0 1600 600])
set(f1,'PaperPositionMode','auto')
f2 = figure(3); %Zeitpunkt 2
set(0,'DefaultAxesFontSize',axissize)
set(f2,'Position',[0 0 1600 600])
set(f2,'PaperPositionMode','auto')
f3 = figure(4); %Zeitpunkt 1
set(0,'DefaultAxesFontSize',axissize)
set(f3,'Position',[0 0 1600 600])
set(f3,'PaperPositionMode','auto')
f4 = figure(5); %Zeitpunkt 2
set(0,'DefaultAxesFontSize',axissize)
set(f4,'Position',[0 0 1600 600])
set(f4,'PaperPositionMode','auto')
subplotleftposition = [0.036 0.1 0.42 0.78];
subplotrightposition = [0.545 0.1 0.42 0.78];


%Cropping parameters
cropx = 382;
cropy = 268;
width = 500;
height = 370;


% Einlesen und croppen der Bilder
p_1=imread('175_p.tif'); %Zeitpunkt 1 (Eingeschaltet)
p_2=imread('175_p2.tif'); %Zeitpunkt 2 (Eingeschaltet)
rf=imread('Referenz_p.tif'); %Ausgeschaltet
rect2 =[cropx  cropy  width  height];
% rect3 = [362.5100    2.5100  523.9800  641.9800]
% with control cell
% rect4 =[382.5100  267.5100  493.9800  369.9800]
% without control cell
p_1s=imcrop(p_1,rect2);
p_2s=imcrop(p_2,rect2);
rfs=imcrop(rf,rect2);


% Zeitpunkt 1
[x1,y1,u1,v1,q1,vld1]=SimplePIV(rfs,p_1s,32,8,16); %Calculation of gradient field

%Plot (Gradient Field 1)
set(0,'CurrentFigure',f1)
s11 =subplot(1,2,1);
set(s11,'Position',subplotleftposition)
quiver(x1,y1,-u1,-v1,2)
view(0,270)
title('Apparent Image Shifts At The Beginning Of The Heating Process (Time 1)', 'Fontsize',titlesize)
xlabel('x [px]', 'Fontsize',axissize)
ylabel('y [px]', 'Fontsize',axissize)

%Plot (Quality Factor 1)
s12 = subplot(1,2,2)
set(s12,'Position',subplotrightposition)
%surfc(x1,y1,q1)
contourf(x1,y1,q1)
c1 = colorbar;
set(get(c1,'ylabel'),'string','Quality Factor [-]','Fontsize',axissize)
view(0,270)
title('Quality Factor Of Automated Correlation Routine (Time 1)', 'Fontsize',titlesize)
xlabel('x [px]', 'Fontsize',axissize)
ylabel('y [px]', 'Fontsize',axissize)



% Zeitpunkt 2
[x2,y2,u2,v2,q2,vld2]=SimplePIV(rfs,p_2s,32,8,16);%Calculation of gradient field

%Plot (Gradient Field 2)
set(0,'CurrentFigure',f2)
s21 = subplot(1,2,1)
set(s21,'Position',subplotleftposition)
quiver(x2,y2,-u2,-v2,2)
view(0,270)
title('Apparent Image Shifts At The Beginning Of The Heating Process (Time 2)', 'Fontsize',titlesize)
xlabel('x [px]', 'Fontsize',axissize)
ylabel('y [px]', 'Fontsize',axissize)

%Plot (Quality Factor 2)
s22 =subplot(1,2,2)
set(s22,'Position',subplotrightposition)
contourf(x2,y2,q2)
c2=colorbar;
set(get(c2,'ylabel'),'string','Quality Factor [-]','Fontsize',axissize)
view(0,270)
title('Quality Factor Of Automated Correlation Routine (Time 2)', 'Fontsize',titlesize)
xlabel('x [px]', 'Fontsize',axissize)
ylabel('y [px]', 'Fontsize',axissize)



% Umrechnung in Gradient des Brechungsindexes
T_0 = 20;
n_0 = 1.48-0.00036*T_0;
L=10*10^(-3);
sigma= 50*10^(-3)/width; % [m/pix]
phi=16*sigma;
 


% Zeitpunkt 1

% Skalierung des Verschiebungsvektors mit Kalibrationsergebnis
gradnx1=-(n_0/L)*tan(2*pi/360)*(u1)/delta(3);
gradny1=-(n_0/L)*tan(2*pi/360)*(v1)/delta(3);

% Integration des Gradientenfeldes des Brechungsindexes
n1=IntegrateDisplacements(gradnx1*phi,gradny1*phi,NaN(height,1),n_0*ones(height,1),NaN(1,width),NaN(1,width));


%Plot (Refraction index 1)
set(0,'CurrentFigure',f3)
s31 = subplot(1,2,1)
set(s31,'Position',subplotleftposition)
contourf(x1,y1,n1)
c3=colorbar;
set(get(c3,'ylabel'),'string','Refraction Index [-]','Fontsize',axissize)
view(0,270)
title('Refractive index (Time 1)', 'Fontsize',titlesize)
xlabel('x [pix]', 'Fontsize',axissize)
ylabel('y [pix]', 'Fontsize',axissize)


% Umrechnung in Temperaturfeld
T1=(1.48-n1)./0.00036;


%Plot (Temperature Field 1)
s32=subplot(1,2,2)
set(s32,'Position',subplotrightposition)
contourf(x1,y1,T1)
c4=colorbar;
set(get(c4,'ylabel'),'string','Temperature [°C]','Fontsize',axissize)
view(0,270)
title('Temperature Distribution (Time 1)', 'Fontsize',titlesize)
xlabel('x [pix]', 'Fontsize',axissize)
ylabel('y [pix]', 'Fontsize',axissize)



% Zeitpunkt 2
% Skalierung des Verschiebungsvektors mit Kalibrationsergebnis
gradnx2=-(n_0/L)*tan(2*pi/360)*(u2)/delta(3);
gradny2=-(n_0/L)*tan(2*pi/360)*(v2)/delta(3);

% Integration des Gradientenfeldes des Brechungsindexes (Grösse des Ausschnittes: x-Richtung 524 pixel, y-Richtung 370 pixel)
n2=IntegrateDisplacements(gradnx2*phi,gradny2*phi,NaN(height,1),n_0*ones(height,1),NaN(1,width),NaN(1,width));

%Plot (Refraction index 2)
set(0,'CurrentFigure',f4)
s41 = subplot(1,2,1)
set(s41,'Position',subplotleftposition)
contourf(x2,y2,n2)
c5=colorbar;
set(get(c5,'ylabel'),'string','Refraction Index [-]','Fontsize',axissize)
view(0,270)
title('Refractive index (Time 2)', 'Fontsize',titlesize)
xlabel('x [pix]', 'Fontsize',axissize)
ylabel('y [pix]', 'Fontsize',axissize)

% Umrechnung in Temperaturfeld
T2=(1.48-n2)./0.00036;


%Plot (Temperature Field 2)
s42= subplot(1,2,2)
set(s42,'Position',subplotrightposition)
view(0,270)
contourf(x2,y2,T2)
c6=colorbar;
set(get(c6,'ylabel'),'string','Temperature [°C]','Fontsize',axissize)
view(0,270)
title('Temperature Distribution (Time 2)', 'Fontsize',titlesize)
xlabel('x [pix]', 'Fontsize',axissize)
ylabel('y [pix]', 'Fontsize',axissize)


%Exporting of plots to LaTeX
print(f1,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/Schlieren Analysis/Report/pics/figure1','-dpng')
print(f2,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/Schlieren Analysis/Report/pics/figure2','-dpng')
print(f3,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/Schlieren Analysis/Report/pics/figure3','-dpng')
print(f4,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/mheusser/Dropbox/experimentelle methoden/Schlieren Analysis/Report/pics/figure4','-dpng')