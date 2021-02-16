%seq=BlueFoxCapture(nframes,exposure_time,gain,triggerMode,whitebalance,xoffset,yoffset,width,height);
clear all
cd('D:\Users\student\Desktop\mexBlueFox\');
n_frames=25;
exposure_time=40000;
seq=BlueFoxCapture(n_frames,exposure_time,0,0,2,424,446,160,80);
filename='433.jpg'; % Write the temperature of the liquid crystal
rgb= seq(1:3,:,:,2);

rgb= permute(rgb,[2,3,1]);
temp(:,:)=rgb(:,:,1);
rgb(:,:,1)=rgb(:,:,3);
rgb(:,:,3)=temp(:,:);

%%% To show picture or not
%imshow(rgb,[0,255]);

avgrgb= uint8(mean(rgb,4));

small=avgrgb;
cd('D:\Users\student\Desktop\mexBlueFox\htcalibrate2011\');
imwrite(small,filename,'jpg');
imshow(small,[0,255]);
title(['Temperature : '  num2str(str2double(regexprep(filename,'.jpg',''))/10) 'degC'])