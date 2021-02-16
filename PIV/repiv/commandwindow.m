for i=1:size(files)
A(:,:,i)=imread(files(i).name);
end
imshow(A(:,:,1),[])
Warning: Image is too big to fit on screen; displaying at 67% 
> In imuitools\private\initSize at 72
  In imshow at 259 
Amin=min(A,[],3);
imshow(A(:,:,1)-Amin,[])
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (A(:,:,1)-Amin,A(:,:,10)-Amin,1,[62,62],[15,15],[0,0],[16,16],[],0);
Undefined function 'PIV_base' for input arguments of type 'uint16'.
 
cd ..
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (A(:,:,1)-Amin,A(:,:,10)-Amin,1,[62,62],[15,15],[0,0],[16,16],[],0);
quiver(xgrid,ygrid,uvecs,vvecs,3)
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (flipud(A(:,:,1)-Amin),flipud(A(:,:,10)-Amin),1,[62,62],[15,15],[0,0],[16,16],[],0);
quiver(xgrid,ygrid,uvecs,vvecs,3)
quiver(xgrid,ygrid,medfilt2(uvecs),medfilt2(vvecs),3)
[N,X]=hist(uvecs(:),-8:0.1:8);
plot(X,N)
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (flipud(conv2(A(:,:,1)-Amin,ones(3,3)),flipud(conv2(A(:,:,10)-Amin,ones(3,3)),1,[62,62],[15,15],[0,0],[16,16],[],0);
 [xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (flipud(conv2(A(:,:,1)-Amin,ones(3,3)),flipud(conv2(A(:,:,10)-Amin,ones(3,3)),1,[62,62],[15,15],[0,0],[16,16],[],0);
                                                                                                                                                                          |
Error: Unbalanced or unexpected parenthesis or bracket.
 
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (flipud(conv2(A(:,:,1)-Amin,ones(3,3))),flipud(conv2(A(:,:,10)-Amin,ones(3,3))),1,[62,62],[15,15],[0,0],[16,16],[],0);
Warning: CONV2 on values of class UINT16 is obsolete.
         Use CONV2(DOUBLE(A),DOUBLE(B)) or CONV2(SINGLE(A),SINGLE(B)) instead. 
> In uint16.conv2 at 11 
Warning: CONV2 on values of class UINT16 is obsolete.
         Use CONV2(DOUBLE(A),DOUBLE(B)) or CONV2(SINGLE(A),SINGLE(B)) instead. 
> In uint16.conv2 at 11 
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (flipud(conv2(double(A(:,:,1)-Amin),ones(3,3))),flipud(conv2(double(A(:,:,10)-Amin),ones(3,3))),1,[62,62],[15,15],[0,0],[16,16],[],0);
[N,X]=hist(uvecs(:),-8:0.1:8);
plot(X,N)
imshow(peaks,[])
close all
imshow(cmaps,[])
Warning: Image is too big to fit on screen; displaying at 50% 
> In imuitools\private\initSize at 72
  In imshow at 259 
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (flipud(A(:,:,1)-Amin),flipud(A(:,:,10)-Amin),0,[62,62],[15,15],[0,0],[16,16],[],0);
[N,X]=hist(uvecs(:),-8:0.1:8);
plot(X,N)
quiver(xgrid,ygrid,medfilt2(uvecs),medfilt2(vvecs),3)
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000])
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000],'r')
Error using stream2 (line 50)
STARTX,STARTY must all be the same size.

Error in streamline (line 63)
      verts = stream2(x,y,u,v,sx,sy,options);
 
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000],[500 500])
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000],[500 500],'r')
quiver(xgrid,ygrid,medfilt2(uvecs),medfilt2(vvecs),3,'r')
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000],[500 500],'r')
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000],[500 500])
streamline(xgrid,ygrid,uvecs,vvecs,[500 750 1000],[500 500 500])
pcolor(uvecs)
close all
pcolor(uvecs)
pcolor(vvecs)
pcolor(curl(uvecs,vvecs))
pcolor(curl(medfilt2(uvecs),medfilt2(vvecs)))
