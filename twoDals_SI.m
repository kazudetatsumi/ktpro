%
%This program decomposes a set of spectrum image where each image was taken
%by a narrow energy slit into the images by the individual spectral components,
%as well as displays each resolved spectral component separately,
%based on the modified alternating least-square algorithm (J.-H. Wang et al,
%Analytica Chimica Acta 476 (203) 93.).
%USAGE:
%1. Set the dimension of the data: ddim
%      ddim[0]= number of slices (images)
%      ddim[1]=number of channels (slit positions) of the data
%2. Set the parameter: numcom=number of components included in the data sets
%3. Prepare the data image set, which should be separately
%   saved as tiff or other supported image format.
%4. Input the data into IDL by "read image" switch.
%5. Set the data name so that ximage[i,*,*] corresponds to one data image set.
%6. Then, compile and run!

yst=1;
xst=1;
ydim=180;
xdim=201;
ddim=[140,xdim*ydim];
numcom=3;

ximage=zeros(201,180,ddim(1));
for i=1:ddim(1), ximage(:,:,i) = imread('NiL_stack.tif', 1310+i);
    xorg=squeeze(ximage(:,:,i));
    %figure(1);imagesc(xorg);
end;

for i=1:ddim(1),
%xorig=squeeze(ximage(xst:xst+xdim-1,yst:yst+ydim-1,i));
xorig=squeeze(ximage(:,:,i));
% figure(2);imagesc(xorig);
end;

xdat=zeros(ddim(1),ddim(2));
for i=1:ddim(1),
xdat(i,:)=squeeze(reshape(ximage(:,:,i),1,1,ddim(2)));
%xdat(i,:)=squeeze(reshape(ximage(xst:xst+xdim-1,yst:yst+ydim-1,i),1,1,ddim(2)));
end;
%for i=1:ddim(1), imagesc(reshape(xdat(i,:),xdim,ydim));figure(3);
%end;

bdat=rand(numcom,ddim(2));
adat=(xdat*bdat')*inv(bdat*bdat');

nn=0;
while  nn < 300;
%PRINT, nn

for i=1:numcom, adat(:,i)=adat(:,i)/norm(adat(:,i));
end;

ind=find(adat < 0);
adat(ind)=0.;

wb=zeros(numcom,numcom);
for i=1:numcom,
%negs=where(bdat[*,i] LT 0,count)
%wb[i,i]=count/N_ELEMENTS(bdat[*,i]);
wb(i,i)=exp(-(nn/100));
end;

bdat=inv(adat' *adat+ wb)*(adat' * xdat+wb*bdat);
negs=find(bdat <0);
bdat(negs)=0.;

wa=zeros(numcom,numcom);
for i=1:numcom,
%egs=where(adat[*,i] LT 0,count)
%wa[i,i]=count/N_ELEMENTS(adat[*,i])
wa(i,i)=exp(-(nn/100));
end;

adat=(xdat*bdat'+adat*wa)*inv(bdat*bdat'+wa);
nn=nn+1;
end;

%for i=0,numcom-1
%WINDOW, i+1,tit='result#'+string(i+1),xs=300, ys=200
%PLOT, adat[*,i]
%end

rimage=zeros(numcom,xdim,ydim);
for i=1:numcom,
rimage(i,:,:)=reshape(bdat(i,:),xdim,ydim);
%WINDOW, numcom+i+1, tit='component#'+string(i+1), xs=xdim,ys=ydim
%TV, rimage[i,*,*]
end;

figure(3);plot(adat, 'DisplayName', 'adat', 'YDataSource', 'adat'); 
for i=1:numcom,aimage=squeeze(rimage(i,:,:));
%bimage=squeeze(rimage(2,:,:));
 figure(4);subplot(1, numcom,i);imagesc(aimage); 
end;
%figure(3);imagesc(bimage); 
%FOR i=0,numcom-1 DO BEGIN
%openw, 7, out_dir+'als_result'+string(i+1)+'.txt'
%printf, 7, format='(e12.4,2x)', adat[*,i]
%close, 7
%ENDFOR

%res=MAKE_ARRAY(DIM=[xdim,ydim],/FLOAT)
%for i=1,numcom,
%res=REFORM(rimage[i,*,*],xdim,ydim)
%WRITE_TIFF, 'Image'+string(i+1)+'.tif', res
%end


