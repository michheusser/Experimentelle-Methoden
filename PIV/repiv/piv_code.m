%This script will import the piv images and do the requested post
%processing. After that it will plot some images.

%% image reading
%one has to first define a flag so the function will not read in the images
%every time the script is eveluated. So before reading the image, define
%flag=0;
if flag==0
    folder=('Set 2');

    filelist=dir([folder '/*.tif']);

    for i=1:size(filelist)
        im_data(:,:,i)=imread([folder '/' filelist(i).name]);
    end
    flag=1;
end
%% determining scale
% from one image one can determine that the box with fluid is 1190 pixels
% long, which is equal to 5cm in nature
x_axis=linspace(0,1280/1190*5,1280);
y_axis=linspace(0,1024/1190*5,1280);
relation=5/1190; %pixel size in cm

%% minimal image for subtraction
Im_min=min(im_data,[],3);

%% using provided PIV_base script doing the base calculation of v field

%piv velocity field without any changes on the images
[xgrid_1,ygrid_1,uvecs_1,vvecs_1,peaks_1,valid_1,cmaps_1] = PIV_base (flipud(im_data(:,:,1)-Im_min),flipud(im_data(:,:,10)-Im_min),1,[62,62],[15,15],[0,0],[16,16],[],0);
quiver(relation*xgrid_1,relation*ygrid_1,uvecs_1,vvecs_1,3);
title('Flow representation without any additional processing')
ylabel('[cm]')
xlabel('[cm]')
figure
[N,X]=hist(uvecs_1(:),-8:0.1:8);
title('Histogram of x vector sizes without any additional processing')
plot(X,N);
xlabel('pixel displacement [-]')
ylabel('Frequency of occurence [-]')

%% exploring possibilities to improve the results
%after these first plots are made, one observes that a lot of pixel locking
%and also quite some outliers are present
%to avoid outliers, a median filter is used when doing the plot with quiver
figure
quiver(relation*xgrid_1,relation*ygrid_1,medfilt2(uvecs_1),medfilt2(vvecs_1),3);
ylabel('[cm]')
xlabel('[cm]')
title('Flow representation with use of median filter to avoid outliers')
figure

%after applying median filter on the data, one can observe a lot of
%outliers beeing eliminated, however this does not affect the pixel locking
%effect
% one possibility to avoid pixel locking, is the another function for
% sub-pixel intepolation, which is included in the PIV_base code
[xgrid_1,ygrid_1,uvecs_1,vvecs_1,peaks_1,valid_1,cmaps_1] = PIV_base (flipud(im_data(:,:,1)-Im_min),flipud(im_data(:,:,10)-Im_min),0,[62,62],[15,15],[0,0],[16,16],[],0);
quiver(relation*xgrid_1,relation*ygrid_1,medfilt2(uvecs_1),medfilt2(vvecs_1),3);
title('Flow representation with use of median filtering and using professors sub pixel intepolation function')
ylabel('[cm]')
xlabel('[cm]')
figure
[N,X]=hist(uvecs_1(:),-8:0.1:8);
title('Histogram of x vector sizes with use of median filtering and using professors sub pixel intepolation function')
plot(X,N);
xlabel('pixel displacement [-]')
ylabel('Frequency of occurence [-]')

%using another interpolation function improves the systematic error
%significalntly, however some pixel locking is still present
% to avoid pixel locking, some image blurring with using of image
% convolution with matrix of 4x4 is used, and both sub pixel evaluations
% are used to evaluate the best results

%using gaussian interpolation
[xgrid_2,ygrid_2,uvecs_2,vvecs_2,peaks_2,valid_2,cmaps_2] = PIV_base (flipud(conv2(double(im_data(:,:,1)-Im_min),ones(4,4))),flipud(conv2(double(im_data(:,:,10)-Im_min),ones(4,4))),1,[62,62],[15,15],[0,0],[16,16],[],0);
figure
quiver(relation*xgrid_2,relation*ygrid_2,medfilt2(uvecs_2),medfilt2(vvecs_2),3);
ylabel('[cm]')
xlabel('[cm]')
title('Flow representation with use of median filtering and and image blurring, gaussing sub pixel interp.')
figure
[N2,X2]=hist(uvecs_2(:),-8:0.1:8);
plot(X2,N2);
xlabel('pixel displacement [-]')
ylabel('Frequency of occurence [-]')
title('Histogram of x vector sizes with use of median filtering and and image blurring, gaussing sub pixel interp.')
%using professors function
[xgrid_3,ygrid_3,uvecs_3,vvecs_3,peaks_3,valid_3,cmaps_3] = PIV_base (flipud(conv2(double(im_data(:,:,1)-Im_min),ones(4))),flipud(conv2(double(im_data(:,:,10)-Im_min),ones(4))),0,[62,62],[15,15],[0,0],[16,16],[],0);

figure
quiver(relation*xgrid_3,relation*ygrid_3,medfilt2(uvecs_3),medfilt2(vvecs_3),3);
ylabel('[cm]')
xlabel('[cm]')
title('Flow representation with use of median filtering and and image blurring, professors sub pixel interp.')
figure
[N,X]=hist(uvecs_3(:),-8:0.1:8);
plot(X,N);
xlabel('pixel displacement [-]')
ylabel('Frequency of occurence [-]')
title('Histogram of x vector sizes with use of median filtering and and image blurring, professors sub pixel interp.')

%% Another plot of with coloured vectors
% this is to better represent color vectors. To get better contrast the
% background is set to black

figure
quiverc(relation*xgrid_3,relation*ygrid_3,medfilt2(uvecs_3),medfilt2(vvecs_3),3);
ylabel('[cm]')
xlabel('[cm]')
title('Flow representation with use of median filtering and and image blurring, professors sub pixel interp.')


%% ploting some flow properties

%since our flow is very laminar, the flow can be assumed to be time
%indepentand, so instead of streamlines one can plot an average image of
%128 images and this would get a good representation of flow as such:

im_avg = 1/128*sum(im_data,3) - double(Im_min);
im_avg_max=max(max(im_avg)');
im_bw=im2bw(im_avg(100:945,71:1260)/im_avg_max,0.05);
figure
imshow(im_bw);
title('Representation of fluid flow as average image of 128 frames')

%another property that is interesting to investigate is flow vorticity,
%which is calculated by curl function
figure
curl_1=((curl(medfilt2(uvecs_3),medfilt2(vvecs_3))));
pcolor(curl_1(7:46,4:76));
colorbar

%plotting the velocity field strain
figure
[gradx, grady]=gradient(medfilt2(uvecs_3)*medfilt2(vvecs_3)');
pcolor(gradx);
title('Representation of fluid strain in x direction')
colorbar
figure
pcolor(grady);
title('Representation of fluid strain in y direction')
colorbar


