 Amean=imread('A.tif');
Bmean=imread('B.tif');
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (Amean,Bmean,1,[62,62],[15,15],[0,0],[16,16],[],0);
xygrid(:,:,2)=ygrid;
xygrid(:,:,1)=xgrid;
uvgrid(:,:,2)=medfilt2(vvecs);
uvgrid(:,:,1)=medfilt2(uvecs);
 [xgrid,ygrid,uvecs2,vvecs2,peaks,valid,cmaps] = PIV_base (Amean,Bmean,1,[62,62],[15,15],uvgrid,xygrid,[],0);