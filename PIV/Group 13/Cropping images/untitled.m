%% Croping images from Run2 and Run3
clc
clear all

run2=imread('Run2_0010.tif');
run3=imread('Run3_0010.tif');


%[run2cr rect2]=imcrop(run2)

rect2=[59  571   42   49];
run2cr=imcrop(run2, rect2);
run3cr=imcrop(run3, rect2);
f1=figure(1);

imshow(run2cr, 'tight')

f2=figure(2);

imshow(run3cr, 'tight')

%% Printing the images

print(f1,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/croprun2','-dpng')
print(f2,'/Network/Servers/mlh34-0.ethz.ch/Volumes/01_MXS_RAID/01_StudentData/sckaspar/Dropbox/Experimentelle Methoden/PIV/Report/pics/croprun3','-dpng')
