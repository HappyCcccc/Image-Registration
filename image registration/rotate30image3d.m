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
B= imread('NM01.jpg');
B = rgb2gray(B);


T = imread('NM01.jpg');
T = rgb2gray(T);


xoffset = 300;
yoffset = 250;
dsize = 50;
tsize = 5;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);


[bm,bn] = size(B);

tx = dsize/2;

target = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
[tm,tn] = size(target);


B = double(B)./256.0;
target = double(target)./256.0;




tic
d = 3;
n = 5000;
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

for i = 1:d
    Width(1,i)=1.0;
    Center(1,i)=0.5;
end

maxWidth = 0.0;
maxIndex = 1;
maxD = 0;

Mn = intmax('int64');
Vn = intmax('int64');
Center(maxIndex,:);


for i = 1:n
    
maxWidth = 0.0;
% split the retrangle
    % decide which axis to split on
for j = 1:d
        if(Width(maxIndex,j)>maxWidth)
            maxWidth = Width(maxIndex,j);
        %    maxD = j;
        end
end

A =[];
for j=1:d
   if(Width(maxIndex,j)>0.9*maxWidth) 
    A = [A,j];
   end
end

maxD = datasample(A,1);


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

GVn = 4*(Vn^(1/2))*(d*log(3/(Vn^(2/d))))^(d/4);

if (i==1||flag)
%if (i==1||V(i)~=V(i-1) || M(i)~= M(i-1))
    for j = 1 : m
    % getrho(j) = area(j,1).^(2/d)/(f(j)-Mn+d*GVn);
   
    % getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+d*GVn);
     %getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+GVn);
     getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+GVn).^(4/d);
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

scatter3(Center(:,1),Center(:,2),Center(:,3),2*n+1,f(:,1),'*');
set(gcf,'Renderer','OpenGL');
%plot3(Center(:,1),Center(:,2),Center(:,3));
%plot(Center(:,1), Center(:,2),'.');
end
b
Mn
toc


function y = funct(c,d)
global bm;
global tm;

s = round((bm-tm)*c(1));
t = round((bm-tm)*c(2));
r = round(c(3)*360);

%{
s = 250;
t = 250;
r = 30;
%}

y = regimage(s,t,r);
end



function y = regimage(s,t,r)

global B;
global target;
global tm;
global tn;

target1 = imrotate(target,30,'bilinear');

[t1m,t1n] = size(target1);

Y = B(s+1:s+tm,t+1:t+tn);
Z = imrotate(Y,r,'bilinear');

[zm,zn] = size(Z);

if(t1m<zm)
    im = t1m;
else im = zm;
end

if(t1n<zn)
    in = t1n;
else in = zn;
end

cost = 0;
co = 0;
for i=1:im
    for j=1:in
        if(target1(i,j)~=0&&Z(i,j)~=0)
            cost = cost+(target1(i,j)-Z(i,j))^2;
            co = co+1;
        end
    end
end

y = cost/co;
end

