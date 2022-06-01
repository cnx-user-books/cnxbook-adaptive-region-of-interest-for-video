%% mathscript codes
%% The part where we determine windows to devide the frame into regions:
%W = frame width
win1 = ceil(.2*W);
win2 = ceil(.4*W);
win3 = ceil(.6*W);
win4 = ceil(.8*W);
win5 = W;

%% The part with image subtraction for motion detection:

%array M = delta values for each pixel, after subtracting previous frame
%V = array of the pixel values for the current frame, after averaging RBG pixel values

M = abs(M);

valMEAN = mean(V);                 %mean of columns of un-subtracted image
valmean = mean(valMEAN);
         
rowmean = mean(M);                %mean of subtracted values
matMEAN = mean(rowmean);


mean1 = mean(rowmean(1,1:win1-1))/matMEAN;       %get mean of each region, normalize by subtracted mean
mean2= mean(rowmean(win1:win2-1))/matMEAN;
mean3 = mean(rowmean(win2:win3-1))/matMEAN;
mean4 = mean(rowmean(win3:win4-1))/matMEAN;
mean5 = mean(rowmean(win4:win5))/matMEAN;
maxmean = max([mean1 mean2 mean3 mean4 mean5]);    %max value of the region means, for normalizing

MOTvals = [mean1 mean2 mean3 mean4 mean5]./maxmean;   %normalize the region means, put into vector

for i=1:5                                  %ignore nan results
   if isnan(MOTvals(i))
      MOTvals(i)=0;
end
end

%% The part where we average and actually crop:

mid = ceil(W/2);                 %define range of motion for RIO center
wideL = ceil(mid-2/3*H);
midL = ceil(wideL/2);

%average results from detection methods
ROIvals = (MOTvals+EDGEvals+FOCUSvals)./3

dev = ceil(2/3*H)

%use region values to weight a deviation from center for ROI
boxmid=sum([ROIvals(1)*-wideL ROIvals(2)*-midL ROIvals(3)*0 ROIvals(4)*midL ROIvals(5)*wideL])*(wideL/(wideL+midL)); 
boxmidpic = ceil(mid + boxmidold);

left = boxmidpic-dev
right = boxmidpic+dev

%make sure the frame is within bounds
if (left <1)
left = 1
right = left + 2*dev
end
if (right > W-1)
right = W-1
left = right - 2*dev
end

%Mout = zeros(size(P,1),size(P,2));

%Mout(:,left-3:left+3) = 10000000000;

%Mout(:,right-3:right+3) = 10000000000;

%P is un-processed image array
Mout = P(:,left:right);  %output cropped image

%% The part with the 2-D convolution for edge detection:
sobel1 = [-1,0,1;-2,0,2;-1,0,1];
sobel2 = [1,2,1;0,0,0;-1,-2,-1];

dims = size(img);
segwidth = ceil(.2*dims(2));

grad1 = conv2(img,sobel1);
grad2 = conv2(img,sobel2);

gradmag = sqrt(grad1.^2 + grad2.^2);
%graddir = atan(grad1/grad2);

edgemask = (gradmag > t);

vecsums = sum(edgemask(3:dims(1)-2,3:dims(2)-2));

fivesums = [sum(vecsums(1:segwidth));sum(vecsums(segwidth+1:2*segwidth));sum(vecsums(2*segwidth+1:3*segwidth));sum(vecsums(3*segwidth+1:4*segwidth));sum(vecsums(4*segwidth+1:dims(2)-4))];

EDGEvals = fivesums/(max(fivesums));

for i=1:5
   if isnan(EDGEvals(i))
       EDGEvals(i)=0;
   end
end

edgemask(:,segwidth) = 1;
edgemask(:,2*segwidth)= 1;
edgemask(:,3*segwidth) = 1;
edgemask(:,4*segwidth) = 1;

%% The part with focus detection:
%function av=scan(lum)

[m n]=size(lum);

diff=abs(m-n);
cut=ceil(.1*m);
width=m-2*cut;
segments=5;
diff=floor(diff/segments);

for i=1:segments
  
    frame=lum(cut:m-cut-1,(i-1)*diff+1:width+(i-1)*diff);

    pwr_spct=abs(fftshift(fft2(frame))).^2;

    
    [u v]=size(pwr_spct);

    [X Y]=meshgrid(-u/2:u/2-1, -v/2:v/2-1);

    [theta, rho,z]=cart2pol(X,Y,pwr_spct);

    rho=round(rho);
  skip=floor((u/2-1)/10);
  rad_av=zeros(1,skip);
  

  for j=0:skip;
      
     rad_av(j+1)=mean(pwr_spct(find(rho==j)));
  end
   
    f=[0:10:u/2-1];
    x=log(f);
    y=log(rad_av);
    x(1)=[];
    y(1)=[];
    vals=polyfit(x,y,1);
    fit(i)=vals(1);
end


raw=fit;

for i=1:segments
    fit(i)=abs(abs(fit(i))-2);
end 


for i=1:segments
    fit(i)=fit(i)-min(fit);
end

fit=fit./max(fit);

fit=abs(fit-1);
%fit=fit./max(fit);

%fit=abs(fit-1);
%fit=fit./max(fit);

r=rem(segments,5);

l=(segments-r)/5;

for i=1:segments
    if isnan(fit(i))
        fit(i)=0;
    end
end

for i=1:5;
   av(i)=mean(fit((i-1)*l+1:i*l));
end
