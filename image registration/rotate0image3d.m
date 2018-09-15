global B ;
global bm;
global bn;
global tm;
global tn;
global target;
B= imread('NM10.jpg');
B = rgb2gray(B);


T = imread('NM10.jpg');
T = rgb2gray(T);


xoffset = 3000;
yoffset = 2500;
dsize = 500;
tsize = 50;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);

tx = dsize/2;
[bm,bn]=size(B);

target = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
%target = T(150:175,150:170);
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
       %     maxD = j;
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
f(maxIndex,1)= qufun(a,d);


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
NF1=qufun(NCE1,d);
NF2 = qufun(NCE2,d);

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

function y = qufun(c,d)
%n = size(c);
sum=0;
for i = 1:d
    sum = sum+ (c(i)).^2;
end
y = sum;
end

function y = rafunct(c,d)
%n = size(c);
nc=zeros(d);
sum=30;
B=4.0;
for i = 1:d
    nc(i)=2*B*c(i)-B-0.1;
    sum = sum+ nc(i).^2-10*cos(2*pi*nc(i));
end
  y=sum;
%y = sum*0.001;
end

%interpolators interp2
function y = funct1(c,d)
global bm;
global tn;


for i = 1
s = round((bm-tn)*c(i));
t = round((bm-tn)*c(i+1));
%r = round(360*c(i+2));
r = 360*c(i+2);
end

%{
s = c(1);
t = c(2);
r = c(3);
%}
%y = s+t;
%{
s = 250;
t = 250;
r = 30;
%}
y = regimage(s,t,r);
end




function y = regimage(s,t,r)
global B;
global tm;
global tn;
global target;


Y = B(s+1:s+tm, t+1:t+tn);
Y = imrotate(Y,r,'bilinear');
%Y = imrotate(Y,r,'crop');
%[ym,yn]=size(Y);

 cost = 0.0;
        for i=1:tm
            for j=1:tn
                if(Y(i,j)~=0)
                cost = cost + (target(i,j)-Y(i,j))^2;
                end
            end;
        end;
        
y = cost;
end




