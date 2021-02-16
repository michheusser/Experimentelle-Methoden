
files=dir('*.tif')


for i=1:size(files)
    A(:,:,i)=imread(files(i).name);
end

Amin=min(A,[],3);

imshow(A(:,:,1)-Amin,[])

[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (A(:,:,1)-Amin,A(:,:,7)-Amin,1,[48,48],[8,8],[0,0],[16,16],[],0);


quiver(xgrid,ygrid,uvecs,vvecs,3)

[N,X]=hist(uvecs(:),-8:0.1:8);
plot(X,N)

imshow(cmaps,[])





quiver(xgrid,ygrid,uvecs,vvecs,3)
streamline(xgrid,ygrid,uvecs,vvecs,[1000 500])
streamline(xgrid,ygrid,uvecs,vvecs,[500 1000])
help streamline
 streamline  Streamlines from 2D or 3D vector data.
    streamline(X,Y,Z,U,V,W,STARTX,STARTY,STARTZ) creates streamlines
    from 3D vector data U,V,W. The arrays X,Y,Z define the coordinates for
    U,V,W and must be monotonic and 3D plaid (as if produced by MESHGRID). 
    STARTX, STARTY, and STARTZ define the starting positions of the stream
    lines.
    
    streamline(U,V,W,STARTX,STARTY,STARTZ) assumes 
          [X Y Z] = meshgrid(1:N, 1:M, 1:P) where [M,N,P]=SIZE(U). 
    
    streamline(XYZ) assumes XYZ is a precomputed cell array of vertex
        arrays (as if produced by STREAM3).
    
    streamline(X,Y,U,V,STARTX,STARTY) creates streamlines from 2D
    vector data U,V. The arrays X,Y define the coordinates for U,V and
    must be monotonic and 2D plaid (as if produced by MESHGRID). STARTX
    and STARTY define the starting positions of the streamlines. A vector
    of line handles is returned.
    
    streamline(U,V,STARTX,STARTY) assumes 
          [X Y] = meshgrid(1:N, 1:M) where [M,N]=SIZE(U). 
    
    streamline(XY) assumes XY is a precomputed cell array of vertex
        arrays (as if produced by STREAM2).
    
    streamline(AX,...) plots into AX instead of GCA.
 
    streamline(...,OPTIONS) specifies the options used in creating
    the streamlines. OPTIONS is specified as a one or two element vector
    containing the step size and maximum number of vertices in a stream
    line.  If OPTIONS is not specified the default step size is 0.1 (one
    tenth of a cell) and the default maximum number of vertices is
    10000. OPTIONS can either be [stepsize] or [stepsize maxverts].
    
    H = streamline(...) returns a vector of line handles.
 
    Example:
       load wind
       [sx sy sz] = meshgrid(80, 20:10:50, 0:5:15);
       h=streamline(x,y,z,u,v,w,sx,sy,sz);
       set(h, 'Color', 'red');
       view(3);
 
    See also stream3, stream2, coneplot, isosurface, smooth3, subvolume,
             reducevolume.

    Reference page in Help browser
       doc streamline

streamline(xgrid,ygrid,uvecs,vvecs,1000,500)
streamline(xgrid,ygrid,uvecs,vvecs,[200:100:1000],[500 500 500 500 500 500 500 500])
Error using stream2 (line 50)
STARTX,STARTY must all be the same size.

Error in streamline (line 63)
      verts = stream2(x,y,u,v,sx,sy,options);
 
streamline(xgrid,ygrid,uvecs,vvecs,[200:100:1000],[500 500 500 500 500 500 500 500 500])
pcolor(uvec)
Undefined function or variable 'uvec'.
 
pcolor(uvecs)
colormap jet
for i=1:size(files)
A(:,:,i)=flipud(imread(files(i).name));
end
Error using imread (line 369)
File "test_0001.tif" does not exist.
 
cd Raw8\
 for i=1:size(files)
A(:,:,i)=flipud(imread(files(i).name));
end
Amin=min(A,[],3);
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (A(:,:,1)-Amin,A(:,:,7)-Amin,0,[64,64],[8,8],[0,0],[16,16],[],0);
quiver(xgrid,ygrid,uvecs,vvecs,3)
Undefined function 'PIV_base' for input arguments of type 'uint16'.
 
cd ..
[xgrid,ygrid,uvecs,vvecs,peaks,valid,cmaps] = PIV_base (A(:,:,1)-Amin,A(:,:,7)-Amin,0,[64,64],[8,8],[0,0],[16,16],[],0);
quiver(xgrid,ygrid,uvecs,vvecs,3)
