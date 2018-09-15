% deal with target image and floating image
%{
choose a piece of image (250:300,250:300)(name it target, but not real target). 
Then rotated target 30 degree into target1. 
Target1 is the real target we want and we are supposed to calculate the cost 
between target1 and B. I choose same size with target in B and name it Y. 
Rotate Y with degree r into Z. Calculate the cost between target1 and Z.

There are black part in target1 and Z.
 How should I deal with the black part and get the exact cost I want.

I  calculate the average cost for all pixes because just pixes 
which are not black counts. For two different rotation, 
the total number of pixes which counts is different. 
But the optimized point I got by the program is 
(0.83,0.35,0.83) which not the (0.55,0.55,0.083). I don?t know the reason yet.
average cost = 0.0027 
%}


global B ;
global bm;
global bn;
global tm;
global tn;
global target;
global Mn;
global optimage;

%{
B= imread('NM10.jpg');
B = rgb2gray(B);


T = imread('NM10.jpg');
T = rgb2gray(T);


xoffset = 300;
yoffset = 250;
dsize = 500;
tsize = 50;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);


[bm,bn] = size(B);

tx = dsize/2;

%jmc
target1 = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
target = imrotate(target1,30,'bilinear');


[tm,tn] = size(target);


B = double(B)./256.0;
target = double(target)./256.0;
%}

%jmc
%B = imread('Brain.png');
%B = imread('NM01.jpg');

A = imread('dog2.png');
A = rgb2gray(A);
figure
imshow(A);
A = double(A)./256.0;
B = imread('dog.jpg');
B = rgb2gray(B);
figure
imshow(B);

B = double(B)./256.0;
[bm,bn] = size(B);

BR = imrotate(B,90,'bilinear');
[bm1,bn1]=size(BR);

target = BR(300:350,300:350);
%target = BR(100:150,150:200);
[tm,tn] = size(target);

Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')
%pause(50)

tic
d = 3;
n = 4500;
Corner = zeros(1,d);% matrix for corner node
Width = zeros(1,d);% matrix for width and length
Center = zeros(1,d);% matrix for center node
area = ones(1,1);% matrix for retrangle area
f = zeros(1,1);% function value for center points
getrho = zeros(1,1);
a = zeros(1,d);
b = zeros(1,d);
V = zeros(1,n);
M = zeros(1,n);

i=0;
%xbest = rand(d,1)
xbest = [0.5,0.5,0.5]
fbest = funct(xbest,d);

for i=1:n
temp = rand(d,1);
fb = funct(temp,d);
if fb < fbest
    fbest = fb;
    xbest = temp;
end
end

x = xbest;
toc
s = x(1);
t = x(2);
r = x(3);
Z = imrotate(B,r*180,'bilinear');
xtrans = s*(bm-tm)
ytrans = t*(bn-tn)
rangle = r*180
cost=funct(x,d)
Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z2 = Z1(1:tm,1:tn);
imshowpair(target, Z2, 'montage')
pause(5)
%Z = imrotate(B,30,'bilinear');
%Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')

%regimage2(299.0/(bm-tm),299.0/(bn-tn),30.0/180.0)
%regimage2(b(1),b(2),b(3))

function y = funct(c,d)
global bm;
global tm;

%jmc

s = c(1);
t = c(2);
r = c(3);

%{
s = 250;
t = 250;
r = 30;
%}

y = regimage(s,t,r);
end


function y = gradient(c)
global bm;
global tm;
global bn;
global tn;

%jmc

s = c(1);
t = c(2);
r = c(3);

%{
s = 250;
t = 250;
r = 30;
%}

y = [(regimage(s+1/(bm-tm),t,r)-regimage(s,t,r))/(1/(bm-tm))
     (regimage(s,t+1/(bn-tn),r)-regimage(s,t,r))/(1/(bn-tn))
     (regimage(s,t,r+1/180)-regimage(s,t,r))/(1/(1/180))];
end

function y = regimage(s,t,r)

global B;
global target;
global tm;
global tn;
%jmc
global bm;
global bn;
global Mn;
global optimage;

Z = imrotate(B,r*180,'bilinear');
%Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z2 = Z1(1:tm,1:tn);


[cost,nn] = sumsqr(target-Z2);
cost = sqrt(cost);

if(cost<Mn)
   optimage = Z2; 
end
%size(target)
%size(Z2)

y = cost;

end
function y = regimage2(s,t,r)

global B;
global target;
global tm;
global tn;
%jmc
global bm;
global bn;
global optimage;

s*(bm-tm)
t*(bn-tn)
r*180


Z = imrotate(B,r*180,'bilinear');
%Z2=Z(1+s*(bm-tm):tm+s*(bm-tm),1+t*(bn-tn):tn+t*(bn-tn));
Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z2 = Z1(1:tm,1:tn);

%{
Z = imrotate(B,30,'bilinear');
Z1 = imtranslate(Z,[-299,-199]);
Z2 = Z1(1:tm,1:tn);
%}


imshowpair(target, optimage, 'montage')

[cost,nn] = sumsqr(target-Z2);
cost = sqrt(cost);

diff = target -Z2;
%diff
%target(1:5,1:5)
%Z2(1:5,1:5)

%size(target)
%size(Z2)

y = cost;

end

