function [x,y,u,v,quality,valid]= SimplePIV (I1,I2,wxy,mxy,sxy)

%
%   uses normalized cross correlation method and Gaussian peak fit
%
%   I1,I2:  images 
%
%   wxy:    interrogation window size
%   mxy:    max search size in x and y 
%   sxy:    window spacing increment
%
%----------------------- global variables -----------------------------------

clear global CSZ
clear global QMAT QQ

global CSZ
CSZ= 2*mxy + 1;                     % size of correlation map, always odd

% matrix for Gaussian peak fit
global QMAT QQ
QMAT= [ [ 1 -1 -1  1  1  1]; ...
        [ 1  0 -1  0  0  1]; ...
        [ 1  1 -1  1 -1  1]; ...
        [ 1 -1  0  1  0  0]; ...
        [ 1  0  0  0  0  0]; ...
        [ 1  1  0  1  0  0]; ...
        [ 1 -1  1  1 -1  1]; ...
        [ 1  0  1  0  0  1]; ...
        [ 1  1  1  1  1  1] ];
QQ= QMAT' * QMAT;

%-----------------------------------------------------------------------------

if size(I1) ~= size(I2)
    disp ('Images should be same size ...');
    return
end

if wxy < 1 
    disp ('Wrong window size wxy ...');
    return
end

if mxy < 1 
    disp ('Search range mxy too small ...');
    return
end

if 2*mxy >= wxy
    disp ('Search range mxy too large ...');
    return
end

if sxy < 1
    disp ('Wrong window spacing sxy ...');
    return
end

[height,width]=size(I1);

% make sure image data are floating point
I1=double(I1);
I2=double(I2);

% compute edge coordinates for PIV grid
left= 1 ;
right= width - wxy + 1;
top= 1;
bot= height - wxy + 1;

% number of PIV grid points
gx= floor((right-left)/sxy) + 1;
gy= floor((bot-top)/sxy) + 1;
gsize= [gy,gx];

wsz= wxy-CSZ+1;

% loop over all subwindows
cols= 0;
for xoff=left:sxy:right
    cols= cols + 1;
    rows= 0;
    for yoff=top:sxy:bot
        rows= rows + 1;
    
        % insert grid coordinated
        x(rows,cols)= xoff + wxy/2 - 0.5;
        y(rows,cols)= yoff + wxy/2 - 0.5;

        % Extract subwindows from the images
        A= I1(yoff+mxy:yoff+wxy-1-mxy,xoff+mxy:xoff+wxy-1-mxy);
        B= I2(yoff:yoff+wxy-1,xoff:xoff+wxy-1);

        try
                        
            valid(rows,cols)= -1;
            
%             cmap= normxcorr2 (A,B);
            
            % compute means, sums-of-squares
            sumA= sum(A(:));
            sumA2= sum(A(:).^2);
            sumB= slidesum(B,wsz,wsz);
            sumB2= slidesum(B.^2,wsz,wsz);

            % compute covariance for windows
            fft_A= fft2(A,wxy,wxy);   % padded for size of subimage B
            fft_B= fft2(B);
            corrAB= real(ifft2(conj(fft_A).*fft_B));
            sumAB= corrAB(1:CSZ,1:CSZ);

            % compute (co-)variances
            norm= wsz.^2;
            varA= sumA2 - sumA.*sumA./norm;
            varB= sumB2 - sumB.*sumB./norm;
            covAB= sumAB - sumA.*sumB./norm;
        
            % compute correlation map
            denom= repmat(varA,CSZ,CSZ) + varB;
            cmap= 2*covAB ./ denom;

            valid(rows,cols)= 0;
            
            % Find value and position (xpeak, ypeak) of the correlation peak
            [peak,index]= max(cmap(:));
            [ypeak,xpeak]= ind2sub (size(cmap),index);
            
%            if abs(xpeak-wxy+mxy)<mxy && abs(ypeak-wxy+mxy)<mxy
            if abs(xpeak-mxy-1)<mxy && abs(ypeak-mxy-1)<mxy

                % compute integer shift relative to window center (= true zero shift)
%                u(rows,cols)= xpeak + mxy - wxy;
%                v(rows,cols)= ypeak + mxy - wxy;
                u(rows,cols)= xpeak - mxy - 1;
                v(rows,cols)= ypeak - mxy - 1;
                valid(rows,cols)= 0;
                quality(rows,cols)= peak;
            
                % try to imporve using Gaussian peak interpolation
                [intpeak,shift,ok]= gausspeak (cmap,[xpeak,ypeak]);
                
                if ok==1 && intpeak/peak > 0.95 && intpeak/peak < 2.0
                    % good data: final interpolated shifts superseedes integer shift estimate
%                     u(rows,cols)= shift(1) + mxy - wxy;
%                     v(rows,cols)= shift(2) + mxy - wxy;
                    u(rows,cols)= shift(1) - mxy - 1;
                    v(rows,cols)= shift(2) - mxy - 1;
                    valid(rows,cols)= 1;
                    quality(rows,cols)= intpeak;
                end

            end
            
        catch

        end

    end
end

return

%----------------------------------------------------------------------------------------

function ss= slidesum(a,n,m)
% sliding window summation; averaging window size is [n,m]
[na,ma]= size(a);
aa= [zeros(1,ma+1); [zeros(na,1),a]];
a1= cumsum(aa,1);
a2 = a1(n+1:na+1,:)-a1(1:na+1-n,:);
a3 = cumsum(a2,2);
ss = a3(:,m+1:ma+1)-a3(:,1:ma+1-m);

%----------------------------------------------------------------------------------------

function [peak,shift,ok]= gausspeak (cmap,start)

% (pseudo-)Gaussian peak fit: p(x,y) = exp(a + b*x + c*y + d*x*x + e*x*y + f*y*y)

% last modified 23-02-03

global QMAT QQ

peak= 0;
r1= 0;
r2= 0;
shift= [0,0];

ind= 1;
for j=-1:1
    for i=-1:1
        rhs(ind)= log(max(cmap(start(2)+j,start(1)+i),1.e-12));
        ind= ind + 1;
    end
end
qr= QMAT' * rhs';
coeffs= QQ \ qr;
mmat= [ [2*coeffs(4) coeffs(5)]; [coeffs(5) 2*coeffs(6)]];
mrhs= [-coeffs(2); -coeffs(3)];
qvec= (mmat \ mrhs)';

if norm(qvec) < 1/sqrt(2)
    ok= 1;
    shift= qvec + start;
    peak= exp(coeffs(1) + coeffs(2)*qvec(1) + coeffs(3)*qvec(2) + ...
              coeffs(4)*qvec(1)*qvec(1) + coeffs(5)*qvec(1)*qvec(2) + coeffs(6)*qvec(2)*qvec(2));
else
    shift= NaN*sizeof(start);
    peak= NaN;
    ok= -1;
end
return
