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

global C;
global B;
global bm;
global bn;
global tm;
global tn;
global target;

C = imread('Brain.png');
%C=imread('rs-trail1.png');
%C=imread('a.png');
%C=imread('bark.000.tiff');
C = double(C)./256.0;
%C =C(1:300,1:300);
%jmc

B = imread('Brain1.png');
%B = imread('rs-trail2.png');
%B = imread('b.png');
%B = imread('bark.030.tiff');
B = double(B)./256.0;
%B = B(1:300,1:300);
[bm,bn] = size(C);

%BR = imrotate(B,30,'bilinear');

%target = BR(100:150,150:200);
%[tm,tn] = size(target);
[tm,tn] = size(target);

%Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')
%pause(50)

tic
d = 3;
n = 3500;
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

%for k = 1:100
%    u = 0.01*rand(d,1);

for i = 1:d
    Width(1,i)=1.0;
    Center(1,i)=0.5;
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

%jmc
GVn = 4*(Vn^(1/2))*(d*log(3/(Vn^(2/d))))^(d/4);
%GVn = 100.0;

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

%scatter3(Center(:,1),Center(:,2),Center(:,3),2*n+1,f(:,1),'*'); %jmc
%set(gcf,'Renderer','OpenGL');
%scatter3(Center(:,1),Center(:,2),Center(:,3),[],f(:,1),'filled','MarkerEdgeColor','k');
%pause(0.000000000000000001);
%plot3(Center(:,1),Center(:,2),Center(:,3));
%plot(Center(:,1), Center(:,2),'.');
end
b
Mn
toc
%{
s = b(1);
t = b(2);
r = b(3);
Z = imrotate(B,r*180,'bilinear');
xtrans = s*(bm-tm)
ytrans = t*(bn-tn)
rangle = r*180
Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z2 = Z1(1:tm,1:tn);
%imshowpair(target, Z2, 'montage')
%pause(5)
%Z = imrotate(B,30,'bilinear');
%Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')
%}
%regimage2(149.0/(bm-tm),99.0/(bn-tn),30.0/180.0)
regimage2(b(1),b(2),b(3))

%end

function y = funct(c,d)

%jmc

s = c(1);
t = c(2);
r = c(3);

y = regimage(s,t,r);
end



function y = regimage(s,t,r)

global C;
global B;
global target;
global tm;
global tn;
%jmc
global bm;
global bn;

Z = imrotate(B,r*180,'bilinear');
[zm,zn]=size(Z);
Z1 = imtranslate(Z,[-s*0.5*zm,-t*0.5*zn]);
%Z2 = Z1(1:bm,1:bn);

[cm,cn]=size(C);
%[cost,nn] = sumsqr(C-Z1);
[cost,nn] = sumsqr(C(1:min(cm,zm),1:min(cn,zn))-Z1(1:min(cm,zm),1:min(cn,zn)));
cost = cost/(min(cm,zm)*min(cn,zn));
cost = sqrt(cost);

%size(target)
%size(Z2)

y = cost;

end
function y = regimage2(s,t,r)
global C;
global B;

global tm;
global tn;
%jmc
global bm;
global bn;

s*bm*0.5
t*bn*0.5
r*180


Z = imrotate(B,r*180,'bilinear');
[zm,zn]=size(Z);
Z1 = imtranslate(Z,[-s*0.5*zm,-t*0.5*zn]);
%Z2 = Z1(1:bm,1:bn);

[cm,cn]=size(C);
%[cost,nn] = sumsqr(C-Z1);
[cost,nn] = sumsqr(C(1:min(cm,zm),1:min(cn,zn))-Z1(1:min(cm,zm),1:min(cn,zn)));
cost = cost/(min(cm,zm)*min(cn,zn));
cost = sqrt(cost);


%{
Z = imrotate(B,30,'bilinear');
Z1 = imtranslate(Z,[-299,-199]);
Z2 = Z1(1:tm,1:tn);
%}


%imshowpair(C, Z1, 'montage')
figure
imshowpair(C, Z1,'diff');

figure
imshowpair(C, Z1, 'montage');




%diff
%target(1:5,1:5)
%Z2(1:5,1:5)

%size(target)?
%size(Z2)

y = cost;

end

