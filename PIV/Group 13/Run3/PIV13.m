clc
clear all
close all

%% read picture
files1 = dir('*.tif');
for i = 1:size(files1)
    A(:,:,i)=imcrop(flipud(imread(files1(i).name)), [38    90    1202    816]);
    % A(:,:,i)=imread(files1(i).name);
end
Amin = min(A,[],3);

framerate = 0.5/(10*(50/1280)); % [1/s] Bilder pro Sekunde
pic = 12.8/framerate;


%% PIV script
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (A(:,:,1)-Amin,A(:,:,10)-Amin,1,[48,48],[8,8],[0,0],[16,16],[],0);


%% plot results

f1 = figure(1);
imshow(cmaps)
title('Output Composite Image Of Correlation','fontsize', 18)

f2= figure(2);
quiver(xgrid,ygrid,uvecs,vvecs,1.5)
axis([0 1202 0 816])
title('Velocity Field','fontsize', 14)
xlabel('Pixel','fontsize', 14)
ylabel('Pixel','fontsize', 14)

f3 = figure(3);
%[sx sy] = meshgrid(900*(0:0.2:1),1400*(0:0.1:1));
%streamline(xgrid,ygrid,uvecs,vvecs,sx, sy)
streamslice(xgrid,ygrid,uvecs,vvecs)
axis([0 1202 0 816])
title('Streamlines','fontsize', 14)
xlabel('Pixel','fontsize', 14)
ylabel('Pixel','fontsize', 14)

f4 = figure(4);
contourf(xgrid,ygrid,divergence(xgrid,ygrid,uvecs,vvecs))
shading flat
caxis([-0.4 0.4])
c1 = colorbar
title('Divergence','fontsize', 14)
xlabel('Pixel','fontsize', 14)
ylabel('Pixel','fontsize', 14)
set(get(c1,'ylabel'),'string','Divergence [m/s^2]','Fontsize',14)


f5 = figure(5);
contourf(xgrid,ygrid,curl(xgrid,ygrid,uvecs,vvecs))
shading flat
caxis([-0.4 0.4])
c2 = colorbar
title('Vorticity','fontsize', 14)
xlabel('Pixel','fontsize', 14)
ylabel('Pixel','fontsize', 14)
set(get(c2,'ylabel'),'string','Curl [m/s^2]','Fontsize',14)

%% Histogramm
f6 = figure(6)
set(f6,'Position', [0 1000 1000 400])
set(f6,'PaperPositionMode','auto')% print() exports the same size as figures are shown

s01=subplot(1,2,1)
set(s01,'Position',[0.045 0.13 0.44 0.78])
[c,h] = hist(mod(uvecs(:),1),20);
bar(h,c,'b')
title('Histogram Of Pixel-Displacement in x-Direction Around Each Pixel ','fontsize', 14)
xlabel('Pixel','fontsize', 14)
ylabel('Quantity','fontsize', 14)


s02 = subplot(1,2,2)
set(s02,'Position',[0.55  0.13 0.44 0.78])
[c,h] = hist(mod(vvecs(:),1),20);
bar(h,c,'b')
title('Histogram Of Pixel-Displacement in y-Direction Around Each Pixel ','fontsize', 14)
xlabel('Pixel','fontsize', 14)
ylabel('Quantity','fontsize', 14)

f7 = figure(7);
scatter(uvecs(:),vvecs(:))
title('Scatter Plot Of Pixel Displacements (Velocities) ','fontsize', 14)
xlabel('u Values (x-Direction)','fontsize', 14)
ylabel('v Values (y-Direction)','fontsize', 14)
grid on


%% Exporting to LaTeX

print(f1,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure1_run3','-dpng')
print(f2,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure2_run3','-dpng')
print(f3,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure3_run3','-dpng')
print(f4,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure4_run3','-dpng')
print(f5,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure5_run3','-dpng')
print(f6,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure6_run3','-dpng')
print(f7,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/figure7_run3','-dpng')