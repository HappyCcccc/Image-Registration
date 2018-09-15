% deal with target image and floating image

global B ;
global bm;
global bn;
global tm;
global tn;
global target;
global b4m;
global b4n;
B= imread('NM10.jpg');
B = rgb2gray(B);


T = imread('NM10.jpg');
T = rgb2gray(T);


xoffset = 3000;
yoffset = 2500;
dsize = 500;
tsize = 50;

B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);
B1 = imrotate(B,30,'bilinear');
B4 = imrotate(B,45,'bilinear');
[b4m,b4n]=size(B4);

[bm1,bn1] = size(B1);
bm1 = round(bm1/2);
bn1 = round(bn1/2);
target =B1(bm1:bm1+tsize,bn1:bn1+tsize);

[bm,bn] = size(B);

tx = dsize/2;


[tm,tn] = size(target);


B = double(B)./256.0;
target = double(target)./256.0;




tic
d = 3;
n = 1000;
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
            maxD = j;
        end
end

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
GVn = d*(Vn*log(1/Vn)).^(4/d);

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
%{
s = round((bm-tm)*c(1));
t = round((bm-tm)*c(2));
r = round(c(3)*360);


s = 250;
t = 250;
r = 30;
%}
s = c(1);
t = c(2);
r = c(3);
y = regimage(s,t,r);
end



function y = regimage(s,t,r)

global B;
global tm;
global tn;
global target;
global b4m;
global b4n;


B2 = imrotate(B,r*360,'bilinear');
%[b2m,b2n]=size(B2);
B3 = imtranslate(B2,[-s*b4m,-t*b4n],'linear');


        cost = 0.0;
        for i=1:tm
            for j=1:tn 
                cost = cost + (target(i,j)-B3(i,j))^2;
            end;
        end;
 
y = cost;


end

