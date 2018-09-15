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

%A = imread('dog2.png');
%A = rgb2gray(A);
%figure
%imshow(A);
%A = double(A)./256.0;
B = imread('dog.jpg');
B = rgb2gray(B);
%figure
%imshow(B);

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
fileID = fopen('costofrdog.txt','w');
for k =1:1

tic
d = 3;
n = 4000;
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
u = 0.001*rand(1,1);
    Corner = Corner + u;
for i = 1:d
    Width(1,i)=1.0;
    Center(1,i)=0.5+u;
end

maxWidth = 0.0;
maxIndex = 1;
maxD = 0;

%jmc Mn = 1.0*intmax('int64');
%jmc Vn = 1.0*intmax('int64');
Mn = realmax;
Vn = 1.0;
Center(maxIndex,:);


for i = 1:n
    
maxWidth = 0.0;
% split the retrangle
    % decide which axis to split on
for j = 1:d
        if(Width(maxIndex,j)>maxWidth)
            maxWidth = Width(maxIndex,j);
            maxD = j;
        end
end
%{
A =[];
for j=1:d
   if(Width(maxIndex,j)>0.9*maxWidth) 
    A = [A,j];
   end
end

maxD = datasample(A,1);
%}

%update corner node and width/length information for new retrangles
%NC1 = zeros(d);

Width(maxIndex, maxD)=1/3*maxWidth;   
Center = Corner + Width/2;

area(maxIndex)=1;
for j=1:d
    area(maxIndex) = area(maxIndex)*Width(maxIndex, j);
end

a = Center(maxIndex,:);
f(maxIndex,1)= funct(a,d);


if(area(maxIndex)<Vn)
    Vn = area(maxIndex);
    flag = 1;
end

if(f(maxIndex,1)<Mn)
    Mn = f(maxIndex,1);
    flag = 1;
    b = a;
end

%{
if (area(maxIndex) ~= Vn || f(maxIndex) ~= Mn)
    flag = 1;
else 
    flag = 0;
end
%}
NC1 = zeros(1,d);
NC2 = zeros(1,d);
NW1 = zeros(1,d); 
NW2 = zeros(1,d);
NCE1 = zeros(1,d);
NCE2 = zeros(1,d);
NA1= ones(1,1);
NA2= ones(1,1);
NF1= zeros(1,1);
NF2 = zeros(1,1);

for j = 1:d
    if (j == maxD)
        NC1(j) = 1/3*maxWidth + Corner(maxIndex,j);
        NC2(j) = 2/3*maxWidth + Corner(maxIndex,j);
        NW1(j) = 1/3*maxWidth;
        NW2(j) = 1/3*maxWidth;
    else 
        NC1(j) = Corner(maxIndex,j);
        NC2(j) =  Corner(maxIndex,j);
        NW1(j) = Width(maxIndex,j);
        NW2(j) =  Width(maxIndex,j);
    end
    NA1 = NA1*NW1(1,j);
    NA2 = NA2*NW2(1,j);
end
NCE1=NC1+NW1/2;
NCE2=NC2+NW2/2;
NF1=funct(NCE1,d);
NF2 = funct(NCE2,d);

if(NF1<Mn)
    Mn = NF1;
    b = NCE1;
    flag = 1;
end
 
if(NF2<Mn)
    Mn = NF2;
    b = NCE2;
    flag = 1;
end

if(NA1<Vn)
    Vn = NA1;
    flag = 1;
end

if(NA2<Vn)
    Vn = NA2;
    flag = 1;
end


%update corner and width information
Corner = [Corner;NC1;NC2]; 
Width = [Width;NW1;NW2];
Center = [Center;NCE1;NCE2];
area = [area;NA1;NA2];
f = [f;NF1;NF2];

[m,d] = size(Corner);

%Vn = min(area);
% Mn = min(f);

% get GVn
%GVn = d*(Vn*log(1/Vn)).^(4/d);
GVn = (Vn * log(n)).^(2/d);     % work well
%jmc
%GVn = 4*(Vn^(1/2))*(d*log(3/(Vn^(2/d))))^(d/4);   work well
%GVn = 100.0;

if (i==1||flag)
%if (i==1||V(i)~=V(i-1) || M(i)~= M(i-1))
    for j = 1 : m
    % getrho(j) = area(j,1).^(2/d)/(f(j)-Mn+d*GVn);
   
    % getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+d*GVn);
    %getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+GVn);
      getrho(j) = area(j,1).^(2/d)*log(1+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+GVn).^(4/d);
     %getrho(j) = area(j,1);
    end
    
        [tmp, ind] = sort(getrho,'descend');
        Center = Center(ind,:);
        Corner = Corner(ind,:);
        Width = Width(ind,:);
        f = f(ind);
        area = area(ind);
        maxIndex = 1;  
else
    maxIndex = maxIndex+1;
end
flag = 0;
Center(maxIndex, :);
%scatter3(Center(:,1),Center(:,2),Center(:,3),'*');

%scatter3(Center(:,1),Center(:,2),Center(:,3),2*n+1,f(:,1),'*'); %jmc
%set(gcf,'Renderer','OpenGL');
scatter3(Center(:,1),Center(:,2),Center(:,3),[],f(:,1),'filled','MarkerEdgeColor','k');
%pause(0.000000000000000001);
%plot3(Center(:,1),Center(:,2),Center(:,3));
%plot(Center(:,1), Center(:,2),'.');
end
toc
b
Mn

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 5000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.001;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = b'; niter = 0; dx = inf;

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = gradient(x);
    gnorm = norm(g);
    % take step:
    xnew = x - alpha*g;
    newalpha = alpha;
    
    while (any(xnew(:)<0|xnew(:)>1))|(funct(xnew,d)>funct(x,d))
       newalpha = 0.5*newalpha;
       xnew = x -newalpha*g;
    end
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % plot current point
    %plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
    % update termination metrics
    niter = niter + 1
    dx = norm(xnew-x);
    x = xnew;
    
end

xopt = x

s = x(1);
t = x(2);
r = x(3);
Z = imrotate(B,r*180,'bilinear');
xtrans = s*(bm-tm)
ytrans = t*(bn-tn)
rangle = r*180
cost=funct(x,d)
fprintf(fileID, 'cost = %f\n',cost);

end;

fclose(fileID);
% Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
% Z2 = Z1(1:tm,1:tn);
% imshowpair(target, Z2, 'montage')
% pause(5)
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
global tn;
global bn;

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

